"""SRIM-2013 "SR Module": tables of stopping powers and ranges.

Re-implementation of ``SRModule.exe`` (SRIM-2013, file version SRIM-2012.01)
which reads an ``SR.IN`` control file and writes

* a full table of electronic and nuclear stopping, projected range and
  longitudinal / lateral straggling between E-Min and E-Max, or
* stopping powers only, for a list of energies (E-Min = E-Max = 0).

The ranges come from the integration routine of the SR module (a PRAL-type
treatment of the projected range and its moments, after Biersack), which
steps the energy on SRIM's table grid in 50 sub-steps per interval.  The
arithmetic follows the Visual Basic program, including its Single precision
variables, so that the output files are reproduced digit for digit.

Usage::

    from lunasrim import sr
    sr.run('SR.IN')                    # writes the output file named in SR.IN
    table = sr.StoppingTable.from_file('SR.IN')
    for row in table.rows: ...
"""
import datetime
import math
import os
import shlex

from . import vbmath
from . import data_files
from .stopping import electronic_stopping, nuclear_stopping, VBOverflow
from .tally import F32
from .trimin import _Reader
from .vbfmt import vb_format, vb_cstr, vb_str

VERSION = 'SRIM-2012.01'

# SRIM's energy grid: mantissas of one decade
MANTISSA = [None, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.25, 2.5, 2.75, 3,
            3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0]

UNITS = {1: ' eV/Angstrom ', 2: ' keV/micron   ', 3: ' MeV/mm       ', 4: ' keV/(ug/cm2) ',
         5: ' MeV/(mg/cm2) ', 6: ' keV/(mg/cm2) ', 7: '  eV/(1E15 atoms/cm2)', 8: ' L.S.S. reduced units'}

CRLF = '\r\n'


# ------------------------------------------------------- VB Variant arithmetic
def _int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def vdiv(a, b):
    """a / b with VB's result type (Integer / Integer gives a Single)"""
    if _int(a) and _int(b):
        return F32(a / b)
    return a / b


def vmul(a, b):
    if _int(a) and _int(b):
        return a * b
    return a * b


# ------------------------------------------------------------------ input
class SrElement(object):
    def __init__(self, z, name, stoich, mass):
        self.z = z
        self.name = name
        self.stoich = stoich
        self.mass = mass


_pow = vbmath.pow


class SrError(ValueError):
    """SRModule.exe stops with a run-time error for this input"""


def _sqrt(x, what):
    if x < 0:
        raise SrError('negative %s variance: SRModule.exe stops here with a run-time error '
                      '(no table is produced for this input)' % what)
    return vbmath.sqrt(x)


class SrInput(object):
    """Content of an SR.IN file"""

    def __init__(self):
        self.output_name = 'SR_OUTPUT.txt'
        self.ion_z = 1
        self.ion_mass = 1.008
        self.phase = 0              # 0 solid, 1 gas
        self.density = 1.0          # g/cm3
        self.correction = 1         # compound (Bragg) correction factor
        self.elements = []
        self.units = 5
        self.emin = 10
        self.emax = 10000
        self.energies = []          # list mode (emin = emax = 0)

    @classmethod
    def read(cls, path):
        with open(path, 'rb') as fh:
            return cls.parse(fh.read().decode('latin-1'))

    @classmethod
    def parse(cls, text):
        r = _Reader(text)
        s = cls()
        r.line()
        r.line()
        s.output_name = r.token(string=True)
        r.line()
        s.ion_z = r.num()
        s.ion_mass = r.num()
        r.line()
        s.phase = r.num()
        s.density = r.num()
        s.correction = r.num()
        r.line()
        n = r.num()
        r.line()
        for i in range(n):
            z = r.num()
            name = r.token(string=True)
            st = r.num()
            m = r.num()
            s.elements.append(SrElement(z, name, st, m))
        r.line()
        s.units = r.num()
        r.line()
        s.emin = r.num()
        s.emax = r.num()
        if s.emin == 0 and s.emax == 0:
            while True:
                try:
                    e = r.num()
                except (ValueError, IndexError):
                    break
                if e == 0:
                    break
                s.energies.append(e)
        return s

    def to_string(self):
        o = ['---Stopping/Range Input Data (Number-format: Period = Decimal Point)',
             '---Output File Name', '"%s"' % self.output_name,
             '---Ion(Z), Ion Mass(u)', '%s   %s' % (self.ion_z, self.ion_mass),
             '---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.',
             '%s    %s    %s' % (self.phase, self.density, self.correction),
             '---Number of Target Elements', ' %d ' % len(self.elements),
             '---Target Elements: (Z), Target name, Stoich, Target Mass(u)']
        for e in self.elements:
            o.append('%s   "%s"   %s   %s' % (e.z, e.name, e.stoich, e.mass))
        o += ['---Output Stopping Units (1-8)', ' %d' % self.units,
              '---Ion Energy : E-Min(keV), E-Max(keV)', ' %s    %s' % (self.emin, self.emax)]
        if self.emin == 0 and self.emax == 0:
            o += ['%s' % e for e in self.energies] + ['0']
        return CRLF.join(o) + CRLF


# ------------------------------------------------------------ element data
_SCOEF_RAW = {}


def _scoef_raw(directory=data_files.DATA_DIR):
    """symbol, name (with their padding) and numbers (as VB Input # reads
    them) of every element in SCOEF03.DAT"""
    if directory in _SCOEF_RAW:
        return _SCOEF_RAW[directory]
    tab = {}
    with open(data_files._find('SCOEF03.DAT', directory), encoding='latin-1') as fh:
        for line in fh:
            t = shlex.split(line, posix=True)
            if len(t) < 95:
                continue
            m = shlex.shlex(line, posix=False)
            m.whitespace_split = True
            raw = list(m)
            sym, name = raw[0][1:-1], raw[1][1:-1]
            nums = []
            for x in t[2:]:
                try:
                    nums.append(int(x))
                except ValueError:
                    nums.append(float(x))
            tab[nums[0]] = (sym, name, nums)
    _SCOEF_RAW[directory] = tab
    return tab


# ------------------------------------------------------------------ table
class StoppingTable(object):
    """Stopping and range table of SRModule.exe"""

    def __init__(self, si, strict=True):
        self.si = si
        # strict: fail where SRModule.exe fails.  Otherwise the integration
        # stops when the accumulated energy is within rounding of E-Max (on
        # SRIM's grid E-Max can be missed by 1e-12 and SRModule.exe then
        # integrates backwards and stops with a run-time error).
        self.strict = strict
        self.rows = []          # (E keV, Se, Sn, range A, long. straggling A, lat. straggling A)
        self._setup()

    @classmethod
    def from_file(cls, path):
        t = cls(SrInput.read(path))
        t.calculate()
        return t

    # ---------------------------------------------------------- 0x40ee50
    def _setup(self):
        si = self.si
        self.emin, self.emax = si.emin, si.emax
        if self.emin > self.emax:
            self.emin, self.emax = self.emax, self.emin
        tot = 0
        for e in si.elements:
            if e.stoich != 0:
                tot = e.stoich + tot
        self.stoich_sum = tot
        per = vdiv(100, tot)
        self.frac = []
        mz = 0
        mm = 0
        for e in si.elements:
            f = vdiv(vmul(e.stoich, per), 100)
            self.frac.append(f)
            mz = mz + vmul(e.z, f)
            mm = mm + vmul(e.mass, f)
        self.mean_z = mz
        self.mean_mass = mm
        rho = si.density
        n = (6.022e+23 * rho) / mm                  # atoms/cm3
        self.atom_density = n
        f = {}
        f[1] = 1e-24 * n
        f[2] = 1e-23 * n
        f[3] = 1e-23 * n
        f[4] = 0.060222 / mm
        f[5] = f[4]
        f[6] = 1000 * f[4]
        f[7] = float(F32(1.0 / 10.0))
        z1, m1 = si.ion_z, si.ion_mass
        f[8] = ((m1 + mm) * 0.1 * _pow(_pow(z1, 0.6667000053372932) + _pow(mz, 0.6667000053372932), 0.5)) / \
               (8.462 * ((mz * z1) * m1))
        self.unit_factor = f

    # ---------------------------------------------------------- 0x404760
    def calculate(self):
        si = self.si
        z1, m1 = si.ion_z, si.ion_mass
        els = si.elements
        nel = len(els)
        frac = self.frac
        # energy grid (keV)
        E = [None]
        i = 1
        while i <= 1000:
            j = (i - 1) % 26 + 1
            k = float(math.floor(F32((i - 1) / 26))) - 3
            E.append(MANTISSA[j] * _pow(10.0, k))
            if E[i] > self.emax:
                break
            i += 1
        n = min(i, 1000)
        self.n = n
        # electronic stopping (internal units), Bragg's rule
        se = [None] + [0] * n
        for i in range(1, n + 1):
            try:
                electronic_stopping(z1, m1, els[0].z, E[i], si.phase, corr=si.correction)
            except VBOverflow as ex:
                if self.strict:
                    raise SrError(str(ex))
                n = i - 1               # tolerant: the table ends below the limit
                self.n = n
                if n < 2:
                    raise SrError(str(ex))
                self.emax = E[n - 1]
                break
        for k, e in enumerate(els):
            for i in range(1, n + 1):
                s = electronic_stopping(z1, m1, e.z, E[i], si.phase, corr=si.correction) * 10
                if s < 0:
                    raise ValueError('negative stopping')
                se[i] = (frac[k] * s) + se[i]
        conv = (0.60223 * si.density) / self.mean_mass
        for i in range(1, n + 1):
            se[i] = conv * se[i]
        self.E = E
        self.se = se
        self._integrate(E, se, conv)

    def _integrate(self, E, se, conv):
        si = self.si
        z1, m1 = si.ion_z, si.ion_mass
        els = si.elements
        frac = self.frac
        n = self.n
        Z1 = float(z1)
        PI = 3.1415926
        N = (si.density * 0.6023) / self.mean_mass
        A = []          # M2/M1
        G = []          # 4 M1 M2 / (M1+M2)^2
        a = []          # screening length
        epsf = []       # reduced energy per eV
        snf = []        # nuclear stopping factor
        sum_a = 0.0
        sum_g = 0.0
        for k, e in enumerate(els):
            ak = vdiv(e.mass, m1)
            A.append(ak)
            sum_a = sum_a + (ak * frac[k])
            g = (4.0 * ak) / _pow(1.0 + ak, 2.0)
            G.append(g)
            sum_g = sum_g + (g * frac[k])
            sc = 0.4685 / (_pow(z1, 0.23) + _pow(e.z, 0.23))
            a.append(sc)
            ef = (sc * ak) / ((((1.0 + ak) * Z1) * e.z) * 14.4)
            epsf.append(ef)
            snf.append(((((PI * sc) * sc) * g) * N) / ef)
        nk = len(els)

        def nuclear(energy, c):
            sn_t = 0.0
            sv = [0.0] * nk
            st = [0.0] * nk
            for k in range(nk):
                eps = energy * epsf[k]
                s = ((_pow(eps, 0.21226) * 0.01321) + eps) + (c * vbmath.sqrt(eps))
                s = (0.5 * vbmath.log((1.1383 * eps) + 1.0)) / s
                s = snf[k] * s
                sv[k] = s
                t = 1.0 / ((6.684 / _pow(eps, 1.0494)) + ((0.197 / _pow(eps, 1.6991)) + 4.0))
                t = ((G[k] * t) * snf[k]) / epsf[k]
                st[k] = t
            return sv, st

        # starting point: 1.1 eV
        i350, i1e0, i198, i214 = 1, 0, 2, 0
        en = MANTISSA[i350 + 1] * _pow(10, i1e0)
        sv, st = nuclear(en, 0.19593)
        sn = 0.0
        stt = 0.0
        for k in range(nk):
            sn = sn + (sv[k] * frac[k])
            stt = stt + (st[k] * frac[k])
        idx = 1
        ek = en / 1000.0
        while E[idx] <= ek:
            idx += 1
        se_e = se[idx]
        stot = sn + se_e
        x38 = sum_a / ((se_e / sn) + 1.0)
        x320 = 1.0 / (1.0 - (stt / ((sum_g * en) * sn)))
        rp = (2.0 * en) / ((sum_a * sn) * ((1.0 / x38) + 1.0))
        r2 = ((1.0 + x38) * rp) * rp
        lat = r2 / ((x320 / x38) + 2.0)
        nsub = 50
        emin_ev = (self.emin * 1000) - 0.2
        emax_ev = self.emax * 1000
        rows = []
        last = None
        cnt = 0
        while True:
            lon = (r2 - lat) - (rp * rp)
            if not (en < emin_ev):
                cnt += 1
                if lon < 0:
                    lon = 0.0
                    if lat < 0:
                        lat = 0.0
                rows.append((en / 1000.0, se_e, sn, rp, _sqrt(lon, 'longitudinal'), _sqrt(lat / 2.0, 'lateral')))
                last = en
            i350 = i198
            i1e0 = i214
            i198 = (i350 + 1) % 26
            if i198 == 0:
                i198 = 26
            if i198 == 1:
                i214 = i214 + 1
            lo = MANTISSA[i350] * _pow(10, i1e0)
            hi = MANTISSA[i198] * _pow(10, i214)
            if emax_ev < hi:
                hi = emax_ev * 1
            de = (hi - lo) / nsub
            if not (de > 0):
                if self.strict:
                    raise SrError('the energy integration does not reach E-Max exactly (accumulated '
                                  'rounding); SRModule.exe stops here with a run-time error')
                break
            for j in range(nsub):
                en = en + de
                sv, st = nuclear(en, 0.19594)
                sn = 0.0
                s310 = 0.0
                sdc = 0.0
                s1b0 = 0.0
                for k in range(nk):
                    sn = sn + (frac[k] * sv[k])
                    s310 = s310 + ((frac[k] * A[k]) * sv[k])
                    sdc = sdc + ((frac[k] * A[k]) * st[k])
                    s1b0 = s1b0 + ((frac[k] * _pow(1.0 - A[k], 2.0)) * st[k])
                ek = en / 1000.0
                while E[idx] < ek:
                    idx += 1
                if idx < 2:
                    idx = 2
                elif idx > n:
                    idx = n
                fr = (ek - E[idx - 1]) / (E[idx] - E[idx - 1])
                se_e = se[idx - 1] + (fr * (se[idx] - se[idx - 1]))
                stot = sn + se_e
                x88 = (2.0 * s310) * en
                xb0 = 2.0 * sdc
                x204 = 4.0 * xb0
                x154 = (4.0 * en) * stot
                x290 = de / en
                x260 = 1.0 - (x204 / x154)
                rp = (((((4.0 * en) * en) - ((x88 + sdc) * rp)) * x290) / (x154 - xb0)) + rp
                r2 = r2 + (((2.0 * rp) * de) / stot)
                x108 = x154 - x204
                x14c = (2.0 * x88) - s1b0
                x178 = ((en * x204) * rp) / stot
                if x260 >= 0.01:
                    lat = (((((r2 - (2.0 * lat)) * x14c) - x178) * x290) / x108) + lat
                else:
                    lat = (0.5 * r2) - ((_pow(en, 2.0) * rp) / (x88 - (0.5 * s1b0)))
            if en >= emax_ev:
                break
        if last is None or en > last:
            lon = (r2 - lat) - (rp * rp)
            if lon < 0:
                lon = 0.0
            rows.append((en / 1000.0, se_e, sn, rp, _sqrt(lon, 'longitudinal'), _sqrt(lat / 2.0, 'lateral')))
        self.rows = rows
        self.conv = conv

    # ------------------------------------------------------------- output
    def lines(self, date=None):
        si = self.si
        sc = _scoef_raw()
        date = date or datetime.datetime.now()
        o = []
        eq = ' ' + '=' * 66
        o.append(eq)
        o.append('              Calculation using SRIM-2006 ')
        o.append('              SRIM version ---> ' + VERSION[:18])
        o.append('              Calc. date   ---> ' + date.strftime('%B %d, %Y '))
        o.append(eq)
        o.append('')
        o.append(' Disk File Name = ' + si.output_name)
        o.append('')
        ion_line = ' Ion = ' + sc[si.ion_z][1] + ' [' + str(si.ion_z) + '] , Mass = ' + \
            vb_cstr(float(F32(si.ion_mass)), True) + ' amu'
        o.append(ion_line)
        o.append('')
        o.append(' Density = ' + vb_format(si.density, ' 0.0000E+00') + ' g/cm3' +
                 vb_format(self.atom_density, ' = 0.0000E+00') + ' atoms/cm3')
        if si.phase == 1:
            o.append(' Target is a GAS ')
        o.append(' ======= Target  Composition ========')
        o.append('    Atom   Atom   Atomic    Mass     ')
        o.append('    Name   Numb   Percent   Percent  ')
        o.append('    ----   ----   -------   -------  ')
        msum = 0
        for k, e in enumerate(si.elements):
            msum = msum + (sc[e.z][2][3] * self.frac[k])
        for k, e in enumerate(si.elements):
            f100 = self.frac[k] * 100
            line = '     ' + (sc[e.z][0] + '  ')[:2]
            line += ('%7s' % str(e.z))[-7:]
            line += ('%10s' % vb_format(f100, '000.00'))[-10:]
            line += ('%10s' % vb_format((sc[e.z][2][3] * f100) / msum, '000.00'))[-10:]
            o.append(line + '   ')
        o.append(' ' + '=' * 36)
        # SRModule.exe uses 1 for a compound correction of 0
        o.append(' Bragg Correction = ' + vb_format((si.correction or 1) - 1, '0.00%'))
        o.append(' Stopping Units = ' + UNITS[si.units])
        o.append(' See bottom of Table for other Stopping units ')
        o.append('')
        o.append(ion_line)
        o.append('')
        o.append('   Ion        dE/dx      dE/dx     Projected  Longitudinal   Lateral')
        o.append('  Energy      Elec.      Nuclear     Range     Straggling   Straggling')
        o.append('-----------  ---------- ---------- ----------  ----------  ----------')
        u = self.unit_factor[si.units] / self.conv

        def rset(s, w=7):
            return ('%*s' % (w, s))[:w] if len(s) > w else '%*s' % (w, s)

        def length(v):
            if v < 10000.0:
                return rset(vb_format(v, '#######')) + ' A   '
            for lim, div, unit in ((1e7, 1e4, ' um  '), (1e10, 1e7, ' mm  '), (1e13, 1e10, ' m   ')):
                if v < lim:
                    return rset(vb_format(v / div, '####.00')) + unit
            return rset(vb_format(v / 1e13, '####.00')) + ' km  '

        for (ek, se_e, sn, rp, lon, lat) in self.rows:
            ev = ek * 1000
            if ev < 1000.0:
                es = rset(vb_cstr(ev)) + ' eV  '
            elif ev < 1000000.0:
                es = rset(vb_format(ev / 1000.0, '####.00')) + ' keV '
            elif ev < 1000000000.0:
                es = rset(vb_format(ev / 1000000.0, '####.00')) + ' MeV '
            else:
                es = rset(vb_format(ev / 1000000000.0, '####.00')) + ' GeV '
            line = es + vb_format(u * se_e, '  0.000E+00') + vb_format(u * sn, '  0.000E+00 ')
            line += length(rp) + length(lon) + length(lat)
            o.append(line)
        o.append('-' * 59)
        o.append(' Multiply Stopping by        for Stopping Units')
        o.append(' -------------------        ------------------')
        f = self.unit_factor
        for k in range(1, 9):
            o.append(vb_format(f[k] / f[si.units], '  0.0000E+00               ') + UNITS[k])
        o.append(eq)
        o.append(' (C) 1984,1989,1992,1998,2006 by J.P. Biersack and J.F. Ziegler')
        return o

    def write(self, path, date=None):
        with open(path, 'wb') as fh:
            fh.write(''.join(l + CRLF for l in self.lines(date)).encode('latin-1'))


class StoppingList(object):
    """Stopping powers for a list of energies (SR.IN with E-Min = E-Max = 0),
    0x410000 of SRModule.exe"""

    def __init__(self, si):
        self.si = si
        self.table = StoppingTable(si)          # composition and unit factors
        self.rows = []

    def calculate(self):
        si = self.si
        z1, m1 = si.ion_z, si.ion_mass
        frac = self.table.frac
        E = list(si.energies)
        se = [0] * len(E)
        sn = [0] * len(E)
        for k, e in enumerate(si.elements):
            for i, en in enumerate(E):
                try:
                    s = electronic_stopping(z1, m1, e.z, en, si.phase, corr=si.correction) * 10
                except VBOverflow as ex:
                    raise SrError(str(ex))
                se[i] = (frac[k] * s) + se[i]
            for i, en in enumerate(E):
                s = nuclear_stopping(z1, m1, e.z, e.mass, en) * 10
                sn[i] = (frac[k] * s) + sn[i]
        corr = si.correction
        if not (corr == 1 or corr == 0):
            for i, en in enumerate(E):
                if en < 10000:
                    x = vbmath.sqrt(vdiv(en, m1) / 25)
                    f = 1 + (corr - 1) / (vbmath.exp((x - 7.0) * 1.48) - (-1.0))
                    se[i] = se[i] * f
        self.rows = list(zip(E, se, sn))

    def lines(self, date=None):
        si = self.si
        sc = _scoef_raw()
        f = self.table.unit_factor[si.units]
        o = ['Ion = ' + sc[si.ion_z][1] + ', Mass =' + vb_str(si.ion_mass)]
        t = 'Tgt = ' + sc[si.elements[0].z][0]
        for e in si.elements[1:]:
            t += ' + ' + sc[e.z][0]
        t += (' (Gas)' if si.phase == 1 else ' (Solid)') + ', Density =' + vb_str(si.density) + ' g/cm3'
        o.append(t)
        o.append('Stopping Units =' + UNITS[si.units])
        o.append('Energy(keV)  S(Elec)    S(Nuc)')
        for en, se, sn in self.rows:
            o.append(vb_format(en, '0.000E+00') + vb_format(f * se, '   0.000E+00') +
                     vb_format(f * sn, '  0.000E+00 '))
        return o

    def write(self, path, date=None):
        with open(path, 'wb') as fh:
            fh.write(''.join(l + CRLF for l in self.lines(date)).encode('latin-1'))


def run(sr_in, directory=None, date=None):
    """Run the SR module on an SR.IN file; the output file (named in SR.IN)
    is written next to it (or into ``directory``).  Returns its path."""
    si = SrInput.read(sr_in)
    out_dir = directory or os.path.dirname(os.path.abspath(sr_in))
    path = os.path.join(out_dir, si.output_name)
    t = StoppingList(si) if (si.emin == 0 and si.emax == 0) else StoppingTable(si)
    t.calculate()
    t.write(path, date)
    return path


def stopping_table(ion_z, ion_mass, elements, density, emin, emax, phase=0, correction=1,
                   units=5):
    """Stopping and range table without an SR.IN file.

    elements: list of (Z, stoichiometry, mass) of the target atoms.
    Returns the list of rows (E keV, S_elec, S_nuc, projected range A,
    longitudinal straggling A, lateral straggling A); the stopping powers
    are in the chosen ``units`` (1-8, see UNITS).
    """
    si = SrInput()
    si.ion_z, si.ion_mass, si.density, si.phase, si.correction = ion_z, ion_mass, density, phase, correction
    si.elements = [SrElement(z, '', st, m) for (z, st, m) in elements]
    si.units, si.emin, si.emax = units, emin, emax
    t = StoppingTable(si, strict=False)
    t.calculate()
    u = t.unit_factor[units] / t.conv
    return [(e, u * se, u * sn, rp, lon, lat) for (e, se, sn, rp, lon, lat) in t.rows]
