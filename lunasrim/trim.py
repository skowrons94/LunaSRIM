"""Monte Carlo transport of ions in matter: a Python port of TRIM (SRIM-2013).

The algorithm, the order of the floating point operations and the use of the
random number generator follow TRIM.exe, so that for the same TRIM.IN (and
the same random seed) the same trajectories are obtained.

Indices: atom index 0 is the ion, 1..n are the target elements; layers are
numbered 0..nlay-1 here (TRIM uses 1..nlay).
"""
import math

from .rng import TrimRandom
from .trim_setup import TrimTarget, PI
from .trimin import f32

_sqrt = math.sqrt
_log = math.log
_exp = math.exp
_cos = math.cos
_sin = math.sin
_atan = math.atan


def _atn(x):
    return _atan(x)


class IonState(object):
    """State of a moving particle (TRIM's global ion variables)."""
    __slots__ = ('e', 'x', 'y', 'z', 'cx', 'cy', 'cz', 'sine', 'layer')


class TrimEngine(object):
    """Runs TRIM ion histories for a :class:`~lunasrim.trimin.TrimInput`."""

    def __init__(self, ti, data=None, log=None):
        self.ti = ti
        self.t = TrimTarget(ti, data)
        self.rng = TrimRandom(ti.seed1, ti.seed2)
        self.log = log
        self.cascade = ti.cascade
        self.kp = 1 if ti.cascade == 1 else 0
        t = self.t
        self.e0 = t.e0
        self.ef = 1                       # lowest energy followed (eV)
        self.xmin = 0.0
        self.xmax = t.xmax
        self.binw = (self.xmax - self.xmin) / 100
        self.nlay = t.nlay
        self.cnt518 = 0
        self.grid_g = t.grid_g
        self.grid_eg = _exp(t.grid_g)
        self.ion_number = 0
        # results of the last ion
        self.fate = None
        # TRIM's KP flag: 1 = Kinchin-Pease, 2 = monolayer collision steps
        self.kpflag = 1 if self.cascade in (1, 4, 6) else (2 if self.cascade == 3 else 0)
        # recoil stack (Single precision arrays in TRIM, index 1..STACK-1)
        n = self.STACK + 2
        self.r_type = [0] * n
        self.r_lay = [0] * n
        self.r_cx = [0.0] * n
        self.r_cy = [0.0] * n
        self.r_cz = [0.0] * n
        self.r_sine = [0.0] * n
        self.r_e = [0.0] * n
        self.r_x = [0.0] * n
        self.r_y = [0.0] * n
        self.r_z = [0.0] * n
        self.z_of = [ti.ion_z] + [e.z for e in ti.elements]
        self.ncx = self.ncy = self.ncz = self.nsine = 0.0

    STACK = 200000

    # ------------------------------------------------------------------ helpers
    def rnd(self, tag):
        v = self.rng.random()
        if self.log is not None:
            self.log.write('RND %r @%s\n' % (v, tag))
        return v

    def find_layer(self, x):
        d = self.t.depth
        n = self.nlay
        for j in range(n):
            if x < d[j]:
                return j
        return n - 1

    def choose_atom(self, layer):
        r = self.rnd('45f817')
        fr = self.t.frac[layer]
        n = self.t.nel
        i = 1
        while i <= n:
            r = r - fr[i]
            if r < 0:
                break
            i += 1
        if i > n:
            i = n
        return i

    # ---------------------------------------------------------- electronic loss
    def eloss(self, p, e, length, layer):
        """returns (e, dE, straggling); e may be raised to 0.001"""
        if e < 0.001:
            e = 0.001
        g = self.grid_g
        k = int(math.floor(_log(e) / g)) + 1
        if k > 999:
            k = 999
        if k < 1:
            k = 1
        eg = self.grid_eg
        a = eg ** (k - 1)
        b = eg * a
        tab = self.t.stop[p][layer]
        s0 = tab[k - 1]
        s = s0 + ((tab[k] - s0) * (e - a)) / (b - a)
        de = length * s
        t = self.t
        w = t.strag_low[layer] * e
        bohr = t.bohr[layer]
        if w > bohr:
            m = t.mass[p]
            q = 1 - 1 / (1 + e / ((m * 1000) * 931494)) ** 2
            w = bohr * ((1 / q - 0.5) / (1 / q - 1))
        r = self.rnd('45df02')
        strag = ((r - 0.5) * w) * _sqrt(length)
        if abs(strag) > de:
            strag = (1 if strag > 0 else (-1 if strag < 0 else 0)) * de
        de = de + strag
        if de > e:
            de = e
        return e, de, strag

    # ---------------------------------------------------------------- free path
    def flight_path(self, st):
        t = self.t
        lay = st.layer
        l0 = t.L0[lay]
        if st.x < l0 and st.e == self.e0:
            self.cnt518 = 2
            return l0
        casc = self.cascade
        if casc == 8 or casc == 3:
            return l0
        n = t.N[lay]
        self.fl_l0 = l0
        eps = t.eps_fac[lay] * st.e
        eeg = _sqrt(t.emin_fac[lay] * eps)
        pmax = t.a_ion[lay] / ((eeg + _sqrt(eeg)) + 0.125 * eeg ** 0.1)
        length = 1 / (((PI * pmax) * pmax) * n)
        xnext = st.x + length * st.cx
        self.fl_xnext = xnext
        self.fl_len = length
        if self.find_layer(xnext) != lay:
            self._layer_boundary(st)
        xn = self.fl_xnext
        if xn < self.xmin or xn > self.xmax:
            self._surface(st)
        if st.x < pmax and st.cx < 0:
            self._surface(st)
        length = self.fl_len
        if casc == 2:
            length = length * self.rnd('44b16d')
            return length if length > l0 else l0
        if self.cnt518 < 1:
            return length
        length = length * self.rnd('44b16d')
        self.cnt518 -= 1
        return length if length > l0 else l0

    def _layer_boundary(self, st):
        t = self.t
        lay = st.layer
        x = st.x
        cx = st.cx
        if cx >= 0:
            d = t.depth[lay] - x
        else:
            d = None
            if lay == 0:
                d = x
            if lay > 0:
                d = x - t.depth[lay - 1]
            if d is None:
                d = 0          # VB: Empty variant
        l0 = self.fl_l0
        if d < t.L0[lay]:
            d = l0
        binw = self.binw
        while True:
            if not (abs(cx * self.fl_len) > d):
                break
            self.fl_len = self.fl_len * self.rnd('44b58e')
            self.fl_xnext = abs(cx * self.fl_len)
            if binw > 3 * l0:
                if self.fl_xnext < binw and x > 10 * l0:
                    m = abs(cx)
                    if m < 0.01:
                        m = 0.01
                    q = d / m
                    self.fl_len = q if q < binw else binw
                    break
            else:
                if self.fl_xnext < l0:
                    self.fl_len = l0
                    break
        self.cnt518 = 2

    def _surface(self, st):
        t = self.t
        x = st.x
        cx = st.cx
        if cx >= 0:
            d = self.xmax - x
            if d < 0:
                d = t.depth[st.layer] - x
        else:
            d = x - self.xmin
        casc = self.cascade
        binw = self.binw
        while True:
            if not (abs(cx * self.fl_len) > d):
                break
            self.fl_len = self.fl_len * self.rnd('44b8ff')
            if casc == 1 or casc == 4 or casc == 6:
                if abs(cx * self.fl_len) < binw:
                    self.fl_len = abs(binw / cx)
                    break
            else:
                if abs(cx * self.fl_len) < self.fl_l0:
                    self.fl_len = self.fl_l0
                    break
        if cx > 0:
            if cx * self.fl_len > d:
                self.fl_len = (d + t.N[st.layer]) / cx
        self.cnt518 = 2

    # -------------------------------------------------------------- scattering
    def magic(self, p, i, e, pimp):
        """returns (psi, phi, recoil_theta, recoil_phi, T)"""
        t = self.t
        eps = e * t.pair_eps[p][i]
        b = pimp / t.pair_a[p][i]
        if eps > 10:
            s2 = 1 / (1 + (1 + b * (1 + b)) * ((2 * eps) * b) ** 2)
            cost = 2 * (1 - s2) - 1
        else:
            r = b
            r1 = _log(abs(b * eps) + 1e-08) * -2.7
            if r1 > b:
                r1 = _log(abs(r1 * eps) + 1e-08) * -2.7
                if r1 > b:
                    r = r1
            it = 1
            while True:
                if not (it < 15):
                    break
                if r > 200:
                    r = 200
                ex1 = 0.18175 * _exp(r * -3.1998)
                ex2 = 0.50986 * _exp(r * -0.94229)
                ex3 = 0.28022 * _exp(r * -0.4029)
                ex4 = 0.02817 * _exp(r * -0.20162)
                v = (((ex1 + ex2) + ex3) + ex4) / r
                v1 = -((((v + 3.1998 * ex1) + 0.94229 * ex2) + 0.4029 * ex3) + 0.20162 * ex4) / r
                fr = ((b * b) / r + (v * r) / eps) - r
                fr1 = ((-b * b) / (r * r) + (v + v1 * r) / eps) - 1
                q = fr / fr1
                r = r - q
                it = it + 1
                if it >= 15:
                    break
                if abs(q / r) < 0.0003:
                    break
            roc = (-2 * (eps - v)) / v1
            sqe = _sqrt(eps)
            cc = (0.011615 + sqe) / (0.0071222 + sqe)
            aa = ((2 * eps) * (1 + 0.99229 / sqe)) * b ** cc
            ff = (_sqrt(1 + aa ** 2) - aa) * ((9.3066 + eps) / (14.813 + eps))
            delta = (ff * ((r - b) * aa)) / (1 + ff)
            c = ((b + delta) + roc) / (r + roc)
            c2 = c ** 2
            s2 = 1 - c2
            cost = 2 * c2 - 1
        if cost > 1:
            cost = 1
        sint = _sqrt(1 - cost * cost)
        den = cost + t.pair_m[p][i]
        if den == 0:
            psi = PI / 2
        else:
            psi = _atn(sint / den)
        if psi < 0:
            psi = psi + PI
        if cost == 1:
            rth = PI / 2
        else:
            rth = _atn(sint / (1 - cost))
        phi = (2 * PI) * self.rnd('468d77')
        rphi = phi - PI
        tt = (t.pair_kin[p][i] * s2) * e
        return psi, phi, rth, rphi, tt

    @staticmethod
    def new_direction(psi, phi, cx, cy, cz, sine):
        cpsi = _cos(psi)
        spsi = _sin(psi)
        cphi = _cos(phi)
        sphi = _sin(phi)
        if sine == 0 or abs(cx) == 1:
            t1 = cpsi * cx
            t2 = -spsi * cphi
            t3 = -spsi * sphi
        else:
            q = spsi / sine
            t1 = cpsi * cx + (spsi * cphi) * sine
            t2 = cpsi * cy - q * (((cy * cx) * cphi) - cz * sphi)
            t3 = cpsi * cz - q * (((cz * cx) * cphi) + cy * sphi)
        nrm = _sqrt((t1 ** 2 + t2 ** 2) + t3 ** 2)
        c1 = t1 / nrm
        return c1, t2 / nrm, t3 / nrm, _sqrt(abs(1 - c1 * c1))

    # --------------------------------------------------------------- ion loop
    def run_ion(self):
        t = self.t
        self.ion_number += 1
        log = self.log
        if log is not None:
            log.write('ION %d\n' % self.ion_number)
        st = IonState()
        st.e = self.e0
        ang = (self.ti.angle * 3.1415926) / 180.0
        st.cx = _cos(ang)
        st.sine = _sin(ang)
        st.cy = st.sine
        st.cz = 0
        st.layer = 0
        st.x = 0
        self.cnt518 = 0
        st.y = ((t.L0[0] * st.cy) * self.rnd('4344a7'))
        st.z = ((t.L0[0] * st.cz) * self.rnd('434541'))
        e0 = self.e0
        ef = self.ef
        while True:
            if log is not None:
                log.write('CYC E=%r X=%r Y=%r Z=%r CX=%r CY=%r CZ=%r L=%r\n' % (
                    st.e, st.x, st.y, st.z, st.cx, st.cy, st.cz, st.layer + 1))
            st.layer = self.find_layer(st.x)
            length = self.flight_path(st)
            st.e, de, strag = self.eloss(0, st.e, length, st.layer)
            st.e = st.e - de
            if st.e > e0:
                st.e = e0
            if st.e < ef:
                st.e = ef
            xold = st.x
            st.x = st.x + length * st.cx
            st.y = st.y + length * st.cy
            st.z = st.z + length * st.cz
            st.layer = self.find_layer(st.x)
            if st.x < 0:
                self.fate = 'B'
                return st
            if st.x > t.depth[-1]:
                self.fate = 'T'
                return st
            atom = self.choose_atom(st.layer)
            pimp = _sqrt(self.rnd('434c1c') / ((PI * length) * t.N[st.layer]))
            psi, phi, rth, rphi, tt = self.magic(0, atom, st.e, pimp)
            st.cx, st.cy, st.cz, st.sine = self.new_direction(psi, phi, st.cx, st.cy, st.cz, st.sine)
            self.ncx, self.ncy, self.ncz, self.nsine = st.cx, st.cy, st.cz, st.sine
            st.e = st.e - tt
            if not (st.x < self.xmin or st.x > self.xmax):
                self.recoil(st, atom, tt, rth, rphi, psi, phi)
            if not (st.e > ef):
                self.fate = 'S'
                return st

    # ----------------------------------------------------------------- recoils
    def recoil(self, st, atom, tt, rth, rphi, psi, phi):
        """TRIMMain.Main.STARTRECOIL: deal with the target atom hit by the ion"""
        if tt < self.ef:
            self.on_phonon(st, atom, tt)
            return
        if self.kpflag == 1:
            self.on_kinchin_pease(st, atom, tt)
            return
        rng = self.rng
        saved = (rng.s1, rng.s2)
        self.full_cascade(st, atom, tt, rth, rphi)
        if self.cascade < 4:
            rng.s1, rng.s2 = saved

    def full_cascade(self, st, atom, tt, rth, rphi):
        t = self.t
        ef = self.ef
        rtyp, rlay_, rcx, rcy, rcz, rsin = self.r_type, self.r_lay, self.r_cx, self.r_cy, self.r_cz, self.r_sine
        re, rx, ry, rz = self.r_e, self.r_x, self.r_y, self.r_z
        l0 = t.L0
        l0pi = t.L0pi
        depth_end = t.depth[-1]
        rtyp[1] = atom
        rlay_[1] = st.layer
        if self.cascade < 6:
            self.ncx, self.ncy, self.ncz, self.nsine = self.new_direction(rth, rphi, st.cx, st.cy, st.cz, st.sine)
        rcx[1] = f32(self.ncx)
        rcy[1] = f32(self.ncy)
        rcz[1] = f32(self.ncz)
        rsin[1] = f32(self.nsine)
        rx[1] = f32(st.x)
        ry[1] = f32(st.y)
        rz[1] = f32(st.z)
        self.displaced = 0          # G4fd674
        self.flag60c = 0            # G4fd60c
        if tt > t.ed[atom]:
            self.displaced = 1
            self.on_first_vacancy(st, atom, tt)
            v = f32(tt - t.lattice[atom])
            re[1] = v if not (v < 0) else 0.0
        k = 1
        state = 'LOOP'
        while True:
            if state == 'LOOP':
                rt = rtyp[k]
                rl = rlay_[k]
                lr = l0[rl]
                rx[k] = f32(rx[k] + lr * rcx[k])
                ry[k] = f32(ry[k] + lr * rcy[k])
                rz[k] = f32(rz[k] + lr * rcz[k])
                if rx[k] > depth_end:
                    self.on_recoil_transmitted(k)
                    state = 'NEXT'
                    continue
                if rx[k] < 0:
                    eperp = (rcx[k] ** 2) * re[k]
                    self.on_recoil_surface(k, rt, eperp)
                    if eperp < t.surface[rt]:
                        rcx[k] = -rcx[k]
                        rx[k] = f32(1e-4)
                    else:
                        re[k] = 0.0
                        state = 'NEXT'
                        continue
                rl = self.find_layer(rx[k])
                e, de, strag = self.eloss(rt, re[k], lr, rl)
                if e != re[k]:
                    re[k] = f32(e)
                if de > re[k]:
                    de = re[k]
                self.on_recoil_move(k, rt, rl, de)
                if re[k] < ef:
                    state = 'STOP'
                    continue
                atom = self.choose_atom(rl)
                pimp = _sqrt(self.rnd('4366cd')) * l0pi[rl]
                zr = self.z_of[rt]
                if zr < 21 and rx[k] < 2 * l0[0]:
                    if zr < 14:
                        pimp = pimp * self.rnd('4368e6')
                    else:
                        b0 = _sqrt(self.rnd('43691d')) * l0pi[rl]
                        a0 = _sqrt(self.rnd('43699e')) * b0
                        pimp = a0 + ((b0 - a0) * (zr - 13)) / 7
                psi, phi, rth, rphi, tt = self.magic(rt, atom, re[k], pimp)
                k1 = k + 1
                re[k1] = f32(re[k] - tt)
                re[k1] = f32(re[k1] - de)
                if not (re[k1] > 0):
                    re[k1] = 0.0
                rx[k1] = rx[k]
                ry[k1] = ry[k]
                rz[k1] = rz[k]
                rlay_[k1] = rl
                rtyp[k1] = rt
                c1, c2, c3, c4 = self.new_direction(psi, phi, rcx[k], rcy[k], rcz[k], rsin[k])
                rcx[k1] = f32(c1)
                rcy[k1] = f32(c2)
                rcz[k1] = f32(c3)
                rsin[k1] = f32(c4)
                rtyp[k] = atom
                re[k] = f32(tt)
                rlay_[k] = rl
                c1, c2, c3, c4 = self.new_direction(rth, rphi, rcx[k], rcy[k], rcz[k], 0)
                rcx[k] = f32(c1)
                rcy[k] = f32(c2)
                rcz[k] = f32(c3)
                rsin[k] = f32(c4)
                self.ncx, self.ncy, self.ncz, self.nsine = c1, c2, c3, c4
                k = k1
                if k > self.STACK - 1:
                    k = self.STACK - 1
                if re[k] >= ef:
                    state = 'LOOP'
                else:
                    state = 'STOP'
            elif state == 'STOP':
                self.on_recoil_stop(k)
                state = 'NEXT'
            else:   # NEXT
                k -= 1
                if k <= 0:
                    return
                self.flag60c = 0
                rt = rtyp[k]
                r = self.rnd('437f62')
                if re[k] < t.ed[rt] and rx[k] > (l0[0] * r) * 9:
                    self.on_recoil_replace(k)
                    state = 'NEXT'
                    continue
                if re[k] >= t.ed[rt]:
                    self.on_vacancy(k)
                    if rx[k] > 2 * l0[0]:
                        re[k] = f32(re[k] - t.lattice[rt])
                else:
                    self.flag60c = 1
                state = 'LOOP' if re[k] > ef else 'STOP'

    # tally hooks (overridden by the statistics collector)
    def on_phonon(self, st, atom, tt):
        pass

    def on_kinchin_pease(self, st, atom, tt):
        pass

    def on_first_vacancy(self, st, atom, tt):
        pass

    def on_recoil_transmitted(self, k):
        pass

    def on_recoil_surface(self, k, rt, eperp):
        pass

    def on_recoil_move(self, k, rt, rl, de):
        pass

    def on_recoil_stop(self, k):
        pass

    def on_recoil_replace(self, k):
        pass

    def on_vacancy(self, k):
        pass
