"""SRIM-2013 electronic stopping powers.

This is a re-implementation of the stopping power routines found in
SRIM-2013's ``SRModule.exe`` / ``TRIM.exe`` (reverse engineered from the
binaries).  All energies are in keV, stopping powers are returned in
eV / (1E15 atoms/cm2) unless stated otherwise.

Structure of the SRIM calculation
---------------------------------

* Protons (Z1 = 1): 10 parameter fit for 2 <= E <= 1000 keV/u, Bethe-Bloch
  with shell, Barkas and Bloch corrections above 1200 keV/u and a linear
  interpolation in between.  Below 2 keV/u the stopping is proportional to
  velocity.
* Helium and Lithium ions use the same fitting function with their own
  coefficient sets (valid for 2 .. 2000 keV/u, interpolation up to 2400 keV/u).
* Heavier ions are scaled from helium stopping with an ion dependent and a
  target dependent polynomial in ln(E) (2 .. 30000 keV/u); above 30 MeV/u the
  stopping follows the Bethe formula with a high energy effective charge.
* Bragg's rule is used for compounds, optionally with a compound correction
  factor which fades out with energy.
"""
import math

from . import data_files

# noble gases -- no phase correction is applied to those
NOBLE_GASES = (2, 10, 18, 36, 54, 86)
# elements which are gaseous in their natural state
GAS_ELEMENTS = (1, 2, 7, 8, 9, 10, 17, 18, 35, 36, 53, 54, 85, 86)
PHASE_FACTOR = 0.92


def _poly(coef, x):
    s = 0.0
    for k, c in enumerate(coef):
        s += c * x ** k
    return s


def _fit(c, e):
    """SRIM 10 parameter stopping fit (c is the 10 element coefficient list)."""
    slow = c[0] * e ** c[1] + c[2] * e ** c[3] + c[8] * e ** 0.5
    shigh = (c[4] / e ** c[5]) * (math.log(c[6] / e + c[7] * e) + c[9] / e ** 0.2)
    return slow * shigh / (slow + shigh)


def bethe(z1, z2, e, data=None, phase=0):
    """High energy stopping (Bethe-Bloch with shell, Barkas and Bloch terms).

    e is the energy in keV/u, phase 0 = solid, 1 = gas target.
    """
    t = data_files.element(z2) if data is None else data[z2]
    beta2 = 1.0 - 1.0 / (1.0 + e / 931494.0) ** 2
    lx = math.log(e / 1000.0)
    shell_d = 0.0 if phase == 1 else _poly(t.shell_d, lx)
    shell_c = _poly(t.shell_c, lx)
    if phase == 1:
        mean_i = t.i_gas[0] if z1 == 1 else t.i_gas[1]
        q = 1.35 / e ** 0.4
    else:
        mean_i = t.i_solid
        q = 1.5 / e ** 0.4 + 45000.0 / (z2 * e ** 1.6)
    barkas = (z1 * (0.001 * e) * q) / (0.001 * e + q)
    y = (z1 / 137.036) ** 2 / beta2
    bloch = -y * (1.20206 - y * (1.042 + -0.8549 * y + 0.343 * y * y))
    l = (math.log(1022000.0 * beta2 / (1.0 - beta2)) - beta2 - shell_d - shell_c
         - math.log(mean_i) + barkas + bloch)
    return (z1 * z1 * z2) / (beta2 * 1961.16) * l


def _light_ion(coef, z1, z2, e, emax, data, phase=0):
    """H, He, Li stopping (e in keV/u)."""
    emin = 2.0
    ehigh = emax * 1.2
    if e > ehigh:
        return bethe(z1, z2, e, data, phase)
    ee = min(max(e, emin), emax)
    s = _fit(coef, ee)
    if e > emax:
        sb = bethe(z1, z2, ehigh, data, phase)
        return s + (sb - s) * (e - emax) / (ehigh - emax)
    if e < emin:
        return s * (e / emin) ** 0.5
    return s


def _log10(x):
    return math.log(x) / math.log(10.0)


def _hi_charge(hi, e):
    """high energy effective charge correction (e in keV/u)."""
    d = _log10(e) - hi[1]
    w = hi[3] if d >= 0 else hi[2]
    x = math.exp(-d / w) + 1.0
    return hi[0] + (1.0 - hi[0]) / x + hi[4] * math.exp(-abs(1.3 + hi[1] - _log10(e)))


def electronic_stopping_element(z1, z2, e_per_u, data=None, phase=0):
    """Electronic stopping of ion z1 in pure element z2 (no phase correction).

    e_per_u in keV/u, result in eV/(1E15 atoms/cm2).
    """
    d = data_files.load_scoef() if data is None else data
    t = d[z2]
    if e_per_u > 1e7:
        return 0.0
    if z1 == 1:
        return _light_ion(t.p_coef, 1, z2, e_per_u, 1000.0, d, phase)
    if z1 == 2:
        return _light_ion(t.he_coef, 2, z2, e_per_u, 2000.0, d, phase)
    if z1 == 3:
        return _light_ion(t.li_coef, 3, z2, e_per_u, 2000.0, d, phase)
    # heavy ions: scaled helium stopping
    emin, emax = 2.0, 30000.0
    ec = min(max(e_per_u, emin), emax)
    lne = math.log(ec)
    ion = d[z1]
    s = (_light_ion(t.he_coef, 2, z2, ec, 2000.0, d, phase) * _poly(ion.ion_poly, lne)
         * z1 * z1 / 4.0 * _poly(t.tgt_poly, lne))
    if e_per_u < emin:
        return s * (e_per_u / emin) ** 0.5
    if e_per_u > emax:
        r_hi = _hi_charge(ion.hi, e_per_u) * bethe(2, z2, e_per_u, d, phase)
        r_lo = _hi_charge(ion.hi, emax) * bethe(2, z2, emax, d, phase)
        return s * r_hi / r_lo
    return s


def phase_factor(z2, phase):
    """Low energy gas/solid phase correction factor of SRIM.

    phase = 0 solid target, 1 gas target.
    """
    if z2 in NOBLE_GASES:
        return 1.0
    natural_gas = 1 if z2 in GAS_ELEMENTS else 0
    if natural_gas == phase:
        return 1.0
    return PHASE_FACTOR if natural_gas else 1.0 / PHASE_FACTOR


def fade_correction(factor, e_per_u, z1=1):
    """The phase correction factor is applied fully at low velocity and
    fades out at high velocity (around 1.2 MeV/u for protons, 400 keV/u for
    heavier ions)."""
    x = 1.48 * (math.sqrt(e_per_u / 25.0) - (7.0 if z1 == 1 else 4.0))
    if x > 700:
        return 1.0
    return 1.0 + (factor - 1.0) / (1.0 + math.exp(x))


def electronic_stopping(z1, m1, z2, e_kev, phase=0, data=None, corr=1.0):
    """Electronic stopping of ion (z1, mass m1 [amu]) with energy e_kev [keV]
    in element z2, including SRIM's gas/solid phase correction (and the
    per-element part of the compound correction ``corr``).

    Returns eV/(1E15 atoms/cm2).
    """
    d = data_files.load_scoef() if data is None else data
    if not m1:
        m1 = d[z1].mai_mass
    e = e_kev / m1
    s = electronic_stopping_element(z1, z2, e, d, phase)
    if corr <= 0:
        corr = 1.0
    return s * fade_correction(phase_factor(z2, phase) * corr, e, z1)


def compound_correction(corr, e_per_u):
    """SRIM compound correction (Bragg rule correction) as applied by the SR
    module: full correction at low velocity, fading out at high velocity."""
    if corr == 1 or corr == 0:
        return 1.0
    return fade_correction(corr, e_per_u, 1)


# --------------------------------------------------------------------------
# nuclear stopping
# --------------------------------------------------------------------------
def reduced_energy(z1, m1, z2, m2, e_kev):
    """ZBL reduced energy (e in keV)."""
    return e_kev * 32.53 * m2 / (z1 * z2 * (m1 + m2) * (z1 ** 0.23 + z2 ** 0.23))


def nuclear_stopping(z1, m1, z2, m2, e_kev, snuc=None):
    """SRIM-2013 nuclear stopping (eV/(1E15 atoms/cm2)).

    Uses the pair specific fit coefficients of SNUC03.dat for reduced
    energies below 30 and the high energy limit ln(eps)/(2 eps) above.
    """
    tab = data_files.load_snuc() if snuc is None else snuc
    eps = reduced_energy(z1, m1, z2, m2, e_kev)
    if eps > 30:
        sn = math.log(eps) / (2.0 * eps)
    else:
        c = tab[(z1, z2)]
        sn = 0.5 * math.log(1.0 + c[0] * eps) / (eps + c[1] * eps ** c[2] + c[3] * eps ** 0.5)
    return sn * 8.462 * z1 * z2 * m1 / ((m1 + m2) * (z1 ** 0.23 + z2 ** 0.23))


def universal_nuclear_stopping(z1, m1, z2, m2, e_kev):
    """ZBL universal nuclear stopping (eV/(1E15 atoms/cm2)), as used by TRIM."""
    eps = reduced_energy(z1, m1, z2, m2, e_kev)
    if eps > 30:
        sn = math.log(eps) / (2.0 * eps)
    else:
        sn = math.log(1.0 + 1.1383 * eps) / (2.0 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * eps ** 0.5))
    return sn * 8.462 * z1 * z2 * m1 / ((m1 + m2) * (z1 ** 0.23 + z2 ** 0.23))


# --------------------------------------------------------------------------
# targets with several elements
# --------------------------------------------------------------------------
def compound_stopping(z1, m1, elements, e_kev, phase=0, corr=1.0, data=None):
    """Stopping in a compound.

    elements : list of (z2, m2, atom_fraction)
    returns (S_e, S_n) in eV/(1E15 atoms/cm2) using Bragg's rule.
    """
    d = data_files.load_scoef() if data is None else data
    if not m1:
        m1 = d[z1].mai_mass
    tot = float(sum(f for _, _, f in elements))
    se = sn = 0.0
    for z2, m2, f in elements:
        se += f / tot * electronic_stopping(z1, m1, z2, e_kev, phase, d, corr)
        sn += f / tot * nuclear_stopping(z1, m1, z2, m2, e_kev)
    if e_kev < 10000:
        se *= compound_correction(corr, e_kev / m1)
    return se, sn
