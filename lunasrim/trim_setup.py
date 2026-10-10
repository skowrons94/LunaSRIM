"""Pre-calculation of the target constants used by TRIM (SRIM-2013).

The variable names follow TRIM's internal arrays.  Index 0 of the "atom"
arrays is the ion, indices 1..n are the target elements (in TRIM.IN order).
"""
import math

from . import data_files
from .stopping import electronic_stopping_element, phase_factor
from .trimin import f32

PI = 3.1415926          # TRIM's value of pi


def vbint(x):
    """VB CInt/CLng (round half to even)"""
    return int(round(x))


class Layer(object):
    pass


class TrimTarget(object):
    """All quantities TRIM derives from TRIM.IN before the first ion."""

    def __init__(self, ti, data=None):
        self.ti = ti
        self.data = data_files.load_scoef() if data is None else data
        d = self.data
        nel = len(ti.elements)
        nlay = len(ti.layers)
        self.nel = nel
        self.nlay = nlay
        # atom properties: index 0 = ion
        self.z = [ti.ion_z] + [e.z for e in ti.elements]
        self.mass = [ti.ion_mass] + [e.mass for e in ti.elements]
        self.ed = [0.0] + list(ti.ed)
        self.lattice = [0.0] + list(ti.lattice)
        self.surface = [0.0] + list(ti.surface)
        self.e0 = ti.energy * 1000            # eV
        self.angle = ti.angle * PI / 180.0
        # layer depths
        self.width = [l.width for l in ti.layers]
        self.depth = []
        acc = None
        for w in self.width:
            acc = w if acc is None else w + acc
            self.depth.append(acc)
        self.xmin = 0
        self.xmax = self.depth[-1]
        self.binwidth = (self.xmax - self.xmin) / 100
        self.inv_binwidth = 1 / self.binwidth
        # target composition
        self.phase = [l.phase for l in ti.layers]
        self.bragg = [l.bragg for l in ti.layers]
        self.stoich_sum = []
        self.frac = []                      # frac[L][i] i = 1..nel (index 0 unused)
        for l in ti.layers:
            s = 0.0
            for x in l.stoich:
                s = s + x
            self.stoich_sum.append(s)
        for j, l in enumerate(ti.layers):
            self.frac.append([0.0] + [f32(x / self.stoich_sum[j]) for x in l.stoich])
        # TRIM treats an element as "gas" in every layer as soon as it is
        # present in one gas layer (array D4fd94c)
        self.elem_phase = [0] * (nel + 1)
        for j, l in enumerate(ti.layers):
            if l.phase == 1:
                for i in range(1, nel + 1):
                    if self.frac[j][i] > 0:
                        self.elem_phase[i] = 1
        self.mean_mass = []
        self.mean_z = []
        for j in range(nlay):
            m = 0
            zz = 0.0
            for i in range(1, nel + 1):
                m = m + self.frac[j][i] * self.mass[i]
                zz = zz + self.frac[j][i] * self.z[i]
            self.mean_mass.append(m)
            self.mean_z.append(vbint(zz))
        # per layer constants
        z1, m1 = self.z[0], self.mass[0]
        # Mass*1000 in TRIM: the ion mass is a VB Single, Single*Integer -> Single
        self.mass1000 = [f32(m1 * 1000)] + [m * 1000 for m in self.mass[1:]]
        self.N = []          # atomic density (1/A^3)
        self.L0 = []         # N^-1/3
        self.L0pi = []
        self.m_ratio = []
        self.a_ion = []      # screening length ion-layer
        self.eps_fac = []    # reduced energy factor
        self.emin_fac = []
        self.bohr = []       # Bohr straggling constant
        self.strag_low = []
        self.c724 = []
        self.c728 = []
        for j, l in enumerate(ti.layers):
            n = l.density * 0.6022 / self.mean_mass[j]
            self.N.append(n)
            self.L0.append(n ** f32(-1.0 / 3.0))
            self.L0pi.append(self.L0[j] / math.sqrt(PI))
            self.m_ratio.append(m1 / self.mean_mass[j])
        for j in range(nlay):
            zt = self.mean_z[j]
            a = (0.5292 * 0.8853) / (zt ** 0.23 + z1 ** 0.23)
            self.a_ion.append(a)
            ef = (a * self.mean_mass[j]) / ((zt * z1) * 14.4 * (m1 + self.mean_mass[j]))
            self.eps_fac.append(ef)
            mr = self.m_ratio[j]
            self.emin_fac.append(((5 * ef) * (1 + mr) ** 2) / (4 * mr))
            self.bohr.append((177 * z1) * math.sqrt(zt * self.N[j]))
            self.strag_low.append(((((3.2 * 12.0) * self.N[j]) / self.mass1000[0]) * (float(zt + z1) ** (8.0 / 3.0))) ** 0.5)
            self.c724.append((float(zt) ** (-7.0 / 3.0)) * 0.01)
            self.c728.append(((float(zt) ** (2.0 / 3.0)) * 0.1334) / math.sqrt(self.mean_mass[j]))
        # pairwise constants: atoms 0..nel (0 = ion) x target elements 1..nel
        n_at = nel + 1
        self.pair_m = [[0.0] * n_at for i in range(n_at)]
        self.pair_kin = [[0.0] * n_at for i in range(n_at)]
        self.pair_a = [[0.0] * n_at for i in range(n_at)]
        self.pair_eps = [[0.0] * n_at for i in range(n_at)]
        for i in range(n_at):
            for k in range(1, n_at):
                m = self.mass[i] / self.mass[k]
                self.pair_m[i][k] = m
                self.pair_kin[i][k] = (4.0 * m) / (1 + m) ** 2
                a = (0.5292 * 0.8853) / (float(self.z[k]) ** 0.23 + float(self.z[i]) ** 0.23)
                self.pair_a[i][k] = a
                self.pair_eps[i][k] = a / ((self.z[i] * self.z[k]) * 14.41 * (1 + m))
        self.build_stopping()

    # ------------------------------------------------------------------ stopping
    def build_stopping(self):
        """Electronic stopping tables S[atom][layer][k] in eV/A on the energy
        grid E_k = exp(G)**(k-1) eV, k = 1..1000, G = ln(E0+1)/999."""
        e0 = self.ti.energy
        g = math.log(e0 * 1000 + 1) / 999.0
        self.grid_g = g
        self.stop = []
        for p in range(self.nel + 1):
            z1, m1 = self.z[p], self.mass[p]
            # energy grid in keV/u
            egrid = [0.0] * 1001
            egrid[1] = 0.001 / m1
            eg = math.exp(g)
            for k in range(2, 1001):
                egrid[k] = egrid[k - 1] * eg
            per_layer = []
            for j in range(self.nlay):
                tab = [0.0] * 1000
                for i in range(1, self.nel + 1):
                    f = self.frac[j][i]
                    if f == 0:
                        continue
                    z2 = self.z[i]
                    for k in range(1000):
                        s = self._se10(z1, z2, egrid[k + 1], self.elem_phase[i], self.bragg[j])
                        tab[k] = tab[k] + (s * f) * self.N[j]
                per_layer.append(tab)
            self.stop.append(per_layer)

    def _se10(self, z1, z2, e_u, phase, corr):
        from .stopping import fade_correction
        s = electronic_stopping_element(z1, z2, e_u, self.data, phase)
        if corr <= 0:
            corr = 1.0
        return (fade_correction(corr * phase_factor(z2, phase), e_u, z1) * s) * 10
