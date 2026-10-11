"""Energy loss and energy straggling of an ion beam crossing a layer.

The mean energy follows the SRIM stopping power (electronic + nuclear,
SRModule.exe's ``S(E)`` with the compound correction):

    dE/dx = -S(E)

The energy spread (variance Omega^2) of the beam grows by the electronic
energy-loss straggling W(E) and is stretched or squeezed by the slope of
the stopping power (an ion slower than the mean loses more energy where S
increases with decreasing E):

    d(Omega^2)/dx = W(E) - 2 S'(E) Omega^2,     S' = dS/dE

Two straggling models are given:

* ``trim``: the model of TRIM.exe (``lunasrim.trim``, electronic energy
  loss of each flight path): W = w^2/12 with w = min(k E, Bohr), where
  Bohr = 177 Z1 (Z2 N)^1/2 eV/A^1/2 (w^2/12 = 4 pi e^4 Z1^2 Z2 N, Bohr's
  formula) with TRIM's relativistic factor, and k E the low energy
  limit, k^2 = 3.2*12 N (Z1+Z2)^(8/3) / (1000 M1), E in eV (it only acts
  below a few hundred eV, so in practice this is Bohr's straggling).  Z2 is
  TRIM's mean target Z (rounded to an integer, as TRIM.exe does).
* ``bohr``: Bohr's formula, W = 4 pi e^4 Z1^2 Z2 N, with the exact
  stoichiometric mean Z2.

The nuclear energy loss is included in the mean energy (as its average),
not in the spread: nuclear collisions give a low energy tail, which a TRIM
calculation (transmitted ions, TRANSMIT.txt) shows in full.
"""
import math

from . import sr
from .stopping import VBOverflow, electronic_stopping, nuclear_stopping

E2 = 14.39964            # e^2 in eV*A
BOHR_K = 4 * math.pi * E2 * E2      # 4 pi e^4, eV^2 A^2
FWHM = 2 * math.sqrt(2 * math.log(2))
AVOGADRO = 0.6022        # 1e24/mol, as in TRIM.exe: N (atoms/A^3) = 0.6022 rho / M
THICKNESS_UNITS = ('A', 'nm', 'um', 'mm', 'ug/cm2', 'mg/cm2', '1e15/cm2', 'atoms/cm2')


class StragglingError(ValueError):
    pass


class _Stopping(object):
    """S(E) in eV/A on a fine logarithmic grid, interpolated log-log"""

    def __init__(self, si, emin, emax, per_decade=200):
        si.emin = si.emax = 0
        n = max(20, int(math.ceil(math.log10(emax / emin) * per_decade)) + 1)
        r = (emax / emin) ** (1.0 / (n - 1))
        energies = [emin * r ** i for i in range(n)]
        # Bragg's rule with the compound correction applied once, as in the
        # SR module's tables and in TRIM (the SR module's list mode,
        # sr.StoppingList, applies it a second time)
        table = sr.StoppingTable(si)                  # composition, unit factors
        rows = []
        for e in energies:
            se = sn = 0.0
            for frac, el in zip(table.frac, si.elements):
                try:
                    se += frac * electronic_stopping(si.ion_z, si.ion_mass, el.z, e, si.phase,
                                                     corr=si.correction) * 10
                except VBOverflow as ex:
                    raise StragglingError(str(ex))
                sn += frac * nuclear_stopping(si.ion_z, si.ion_mass, el.z, el.mass, e) * 10
            rows.append((e, se, sn))
        f = table.unit_factor[1]                      # eV/A
        self.table = table
        self.lnE = [math.log(e) for e, se, sn in rows]
        self.lnSe = [math.log(max(f * se, 1e-300)) for e, se, sn in rows]
        self.lnSn = [math.log(max(f * sn, 1e-300)) for e, se, sn in rows]
        self.h = math.log(r)

    def _interp(self, tab, e):
        x = (math.log(e) - self.lnE[0]) / self.h
        i = min(max(int(x), 0), len(tab) - 2)
        t = x - i
        return math.exp(tab[i] + t * (tab[i + 1] - tab[i]))

    def __call__(self, e):
        """(S_elec, S_nuc) in eV/A at e keV"""
        return self._interp(self.lnSe, e), self._interp(self.lnSn, e)

    def total(self, e):
        se, sn = self(e)
        return se + sn

    def slope(self, e):
        """dS/dE (1/A, both in eV)"""
        d = 1e-3
        return (self.total(e * (1 + d)) - self.total(e * (1 - d))) / (2 * d * e * 1000.0)


def to_angstrom(value, unit, density, mean_mass):
    """thickness in A from ``value`` in ``unit`` (see THICKNESS_UNITS)"""
    n = AVOGADRO * density / mean_mass           # atoms/A^3
    if unit == 'A':
        return value
    if unit == 'nm':
        return value * 10.0
    if unit == 'um':
        return value * 1e4
    if unit == 'mm':
        return value * 1e7
    if unit == 'ug/cm2':
        return value * 1e-6 / density * 1e8
    if unit == 'mg/cm2':
        return value * 1e-3 / density * 1e8
    if unit == '1e15/cm2':
        return value * 1e15 / (n * 1e16)
    if unit == 'atoms/cm2':
        return value / (n * 1e16)
    raise StragglingError('unknown thickness unit %r' % unit)


def energy_loss(ion_z, ion_mass, energy, elements, density, thickness, unit='A', phase=0,
                correction=1.0, steps=500):
    """Energy loss and straggling of ions of ``energy`` keV crossing a layer.

    elements: list of (Z, stoichiometry, mass) of the target atoms;
    thickness in ``unit`` (A, nm, um, mm, ug/cm2, mg/cm2, 1e15/cm2,
    atoms/cm2).  Returns a dict with the energy at the exit (keV), the
    energy loss, the straggling (one standard deviation, keV, both models),
    the FWHM, the areal density, and the depth profile ``depth`` (A),
    ``energy``, ``sigma``, ``sigma_bohr`` (keV).  If the ions stop in the
    layer, ``stopped`` is True and the profile ends where the mean energy
    reaches about 1% of the initial energy.
    """
    if not energy > 0:
        raise StragglingError('the ion energy must be positive')
    if not density > 0:
        raise StragglingError('the density must be positive')
    si = sr.SrInput()
    si.ion_z, si.ion_mass = ion_z, ion_mass
    si.density, si.phase, si.correction = density, phase, correction
    si.elements = [sr.SrElement(z, '', st, m) for (z, st, m) in elements if st > 0]
    if not si.elements:
        raise StragglingError('the target has no elements')
    si.units = 1
    stop = _Stopping(si, max(energy * 1e-3, 1e-3), energy * 1.02)
    tab = stop.table
    mean_mass = tab.mean_mass
    x_total = to_angstrom(float(thickness), unit, density, mean_mass)
    if not x_total > 0:
        raise StragglingError('the thickness must be positive')
    # target constants (as TRIM.exe: integer mean Z, N in atoms/A^3)
    tot = sum(st for (z, st, m) in elements if st > 0)
    zz = 0
    for (z, st, m) in elements:
        if st > 0:
            zz = int(round(st / tot * z + zz))
    z_exact = sum(st / tot * z for (z, st, m) in elements if st > 0)
    n = AVOGADRO * density / mean_mass
    bohr = 177 * ion_z * math.sqrt(zz * n)
    k_low = math.sqrt(3.2 * 12.0 * n / (ion_mass * 1000) * (zz + ion_z) ** (8.0 / 3.0))
    w_bohr = BOHR_K * ion_z ** 2 * z_exact * n            # eV^2/A

    def w_trim(e):
        e = e * 1000.0                  # TRIM.exe works in eV
        w = k_low * e
        if w > bohr:
            q = 1 - 1 / (1 + e / (ion_mass * 1000 * 931494)) ** 2
            w = bohr * ((1 / q - 0.5) / (1 / q - 1))
        return w * w / 12.0

    def deriv(e, v, vb):
        s = stop.total(e)
        sp = stop.slope(e)
        return -s / 1000.0, w_trim(e) - 2 * sp * v, w_bohr - 2 * sp * vb   # keV/A, eV^2/A

    e_stop = energy * 0.01
    h = x_total / steps
    x, e, v, vb = 0.0, float(energy), 0.0, 0.0
    se0, sn0 = stop(e)
    prof = {'depth': [0.0], 'energy': [e], 'sigma': [0.0], 'sigma_bohr': [0.0],
            'se': [se0], 'sn': [sn0]}
    eloss_e = eloss_n = 0.0
    stopped = False
    for i in range(steps):
        # RK4 for (E, variance) over one step
        k1 = deriv(e, v, vb)
        k2 = deriv(max(e + h / 2 * k1[0], e_stop / 2), v + h / 2 * k1[1], vb + h / 2 * k1[2])
        k3 = deriv(max(e + h / 2 * k2[0], e_stop / 2), v + h / 2 * k2[1], vb + h / 2 * k2[2])
        k4 = deriv(max(e + h * k3[0], e_stop / 2), v + h * k3[1], vb + h * k3[2])
        de = h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        e_new = e + de
        se, sn = stop(max(e + de / 2, e_stop / 2))
        if e_new <= e_stop:
            stopped = True
            break
        eloss_e += -de * se / (se + sn)
        eloss_n += -de * sn / (se + sn)
        v = max(v + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), 0.0)
        vb = max(vb + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), 0.0)
        e = e_new
        x += h
        se, sn = stop(e)
        prof['depth'].append(x)
        prof['energy'].append(e)
        prof['sigma'].append(math.sqrt(v) / 1000.0)
        prof['sigma_bohr'].append(math.sqrt(vb) / 1000.0)
        prof['se'].append(se)
        prof['sn'].append(sn)
    sigma = math.sqrt(v) / 1000.0
    return {
        'energy_in': float(energy), 'energy_out': e, 'energy_loss': energy - e,
        'electronic_loss': eloss_e, 'nuclear_loss': eloss_n,
        'sigma': sigma, 'fwhm': FWHM * sigma,
        'sigma_bohr': math.sqrt(vb) / 1000.0, 'fwhm_bohr': FWHM * math.sqrt(vb) / 1000.0,
        'thickness': x_total, 'depth_reached': x, 'stopped': stopped,
        'areal_density': n * 1e16 * x_total,             # atoms/cm2
        'mass_thickness': density * x_total * 1e-8 * 1e6,   # ug/cm2
        'mean_energy': _mean(prof['depth'], prof['energy']),
        'atom_density': n * 1e24, 'mean_z': z_exact, 'trim_z': zz,
        'profile': prof,
    }


def _mean(x, y):
    if len(x) < 2 or x[-1] <= 0:
        return y[0]
    s = 0.0
    for i in range(1, len(x)):
        s += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2
    return s / x[-1]
