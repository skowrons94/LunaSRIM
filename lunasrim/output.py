"""Output files of TRIM (SRIM-2013 formats).

The text written here reproduces the files TRIM.exe produces in its
``SRIM Outputs`` directory.  Numbers are formatted with the Visual Basic
conventions (see :mod:`lunasrim.vbfmt`).
"""
import datetime
import math
import os

from . import data_files
from .trimin import f32
from .vbfmt import vb_format, vb_str, vb_cstr, vb_print_num
from .tally import F32

CRLF = '\r\n'

LINE78 = '=' * 78
SRIM_LINE = '============================== SRIM-2013.00 =============================='


def symbol2(z):
    """element symbol right aligned in 2 characters (TRIM's symbol table)"""
    return '%2s' % data_files.element(z).symbol


def energy_label(e_ev):
    """0x45ebf0"""
    if e_ev >= 1e9:
        return vb_cstr(e_ev / 1e9)[:4] + ' GeV'
    if e_ev >= 1e6:
        return vb_cstr(e_ev / 1e6)[:4] + ' MeV'
    if e_ev >= 1e3:
        return vb_cstr(e_ev / 1e3)[:4] + ' keV'
    return vb_cstr(int(round(e_ev)))[:4] + ' eV '


def depth_label(x):
    """0x45eeb0 (x in Angstrom)"""
    sign = ' ' if x >= 0 else '-'
    a = abs(x)
    for lim, unit in ((1e13, ' km'), (1e10, ' m'), (1e7, ' mm'), (1e4, ' um')):
        if a >= lim:
            return sign + vb_str(a / lim)[:5][-4:] + unit
    return sign + vb_str(int(round(a)))[:5][-4:] + ' A'


def cos_string(c, low=-0.9999999):
    """TRIM writes direction cosines through a string, clamped to 0.9999999"""
    s = vb_cstr(c)
    v = float(s) if s not in ('', '-') else 0.0
    if v > 0.9999999:
        return 0.9999999
    if v < low:
        return -0.9999999
    return v


def _rset(s, w):
    """VB RSet: right-align in w characters, keep the leftmost w if longer"""
    return s[:w] if len(s) > w else s.rjust(w)


class TrimOutput(object):
    """Writes the SRIM Outputs files during and after a TRIM run."""

    def __init__(self, engine, directory, date=None):
        self.eng = engine
        self.ti = engine.ti
        self.dir = directory
        os.makedirs(directory, exist_ok=True)
        self.date = date or datetime.datetime.now()
        self.files = {}
        ti = self.ti
        t = engine.t
        names = '+'.join(l.name for l in ti.layers)
        title = ('====== TRIM Calc.= ' + symbol2(ti.ion_z) + '(' + energy_label(ti.energy * 1000) +
                 ') ==> ' + names[:24] + '(' + depth_label(t.depth[-1]) + ') ')
        self.title = (title + '=' * 78)[:78]
        self.header_lines = [
            LINE78,
            '===== TRIMOUT.txt: File of Transmitted / Backscattered / Sputtered Atoms =====',
            '=  This file tabulates the kinetics of ions or atoms leaving the target.     =',
            '=  Column #1: S= Sputtered Atom, B= Backscattered Ion, T= Transmitted Ion.   =',
            '=  Col.#2: Ion Number, Col.#3: Z of atom leaving, Col.#4: Atom energy (eV).  =',
            '=  Col.#5-7: Last location:  X= Depth into target, Y,Z= Transverse axes.     =',
            '=  Col.#8-10: Cosines of final trajectory.                                   =',
            '= *** This data file is in the same format as TRIM.DAT (see manual for uses).=',
            self.title,
            ' Ion  Atom   Energy        Depth       Lateral-Position        Atom Direction      ',
            ' Numb Numb    (eV)          X(A)        Y(A)       Z(A)      Cos(X)  Cos(Y) Cos(Z) ',
        ]
        self.g_b2c = SRIM_LINE
        self.g_b30 = LINE78

    # ----------------------------------------------------------------- files
    def _open(self, key, name):
        f = open(os.path.join(self.dir, name), 'w', encoding='latin-1', newline='')
        self.files[key] = f
        return f

    def _w(self, key, line=''):
        self.files[key].write(line + CRLF)

    def close(self):
        for f in self.files.values():
            f.close()
        self.files = {}

    def start(self):
        """open the per-ion files (TRIMOUT, TRANSMIT, BACKSCAT, SPUTTER, RANGE_3D)"""
        ti = self.ti
        # TRIMOUT.txt (0x43bbc0)
        self._open('out', 'TRIMOUT.txt')
        self._w('out', self.g_b2c)
        self._w('out', self.g_b30)
        for l in self.header_lines[1:]:
            self._w('out', l)
        if ti.file_backscatt == 1:
            self._open('back', 'BACKSCAT.txt')
            self._w('back', self.g_b2c)
            self._w('back', self.g_b30)
            self._w('back', '========== SRIM Outputs\\BACKSCATSRIM Outputs\\ : File of Backscattered Ions =============')
            for l in self.header_lines[2:]:
                self._w('back', l)
        if ti.file_transmit == 1:
            self._open('trans', 'TRANSMIT.txt')
            self._w('trans', self.g_b2c)
            self._w('trans', self.g_b30)
            self._w('trans', '=================  TRANSMIT.txt : File of Transmitted Ions  ==================')
            for l in self.header_lines[2:]:
                self._w('trans', l)
        if ti.file_sputtered == 1:
            self._sputter_header()
        if ti.file_ranges == 1:
            self._range3d_header()

    def _sputter_header(self):
        ti = self.ti
        t = self.eng.t
        self._open('sput', 'SPUTTER.txt')
        w = lambda s='': self._w('sput', s)
        w(self.g_b2c)
        w(self.g_b30)
        w('-' * 78)
        w(' =====================================================================')
        w('                       SPUTTERING   CALCULATION                       ')
        w(' =====================================================================')
        if ti.cascade == 3:
            w('    >> Review plot of SPUTTERING YIELD vs. Surface Binding Energy <<  ')
            w('    >> for how yield changes with SBE (last plot of Menu of Plots)<<  ')
        w(' =====================================================================')
        for l in (' =   AXIS DEFINITIONS: X=Depth, Y,Z= Lateral plane of target surface.=',
                  ' =   (If beam enters target at an angle, this tilt is in Y direction)=',
                  ' =   Shown are: Ion Number                                           =',
                  ' =              Sputtered atom type (atomic name) and energy (eV)    =',
                  ' =                 (Note: SBE is NOT subtracted from atom energy)    =',
                  ' =              Point of ejection relative to ion entrance point.    =',
                  ' =              Trajectional cosines of ejected atom.                =',
                  ' =                  Cosine (X) is negative indicating motion away    =',
                  ' =                     from the target surface.                      =',
                  ' =                  Cosine (Y) and Cosine (Z) are transverse motion. =',
                  ' =====================================================================',
                  ' ==================  CALCULATION   DATA  ============================='):
            w(l)
        d = self.date
        w('  Calculation Date= %02d-%02d-%04d, Time= %02d:%02d:%02d' % (d.month, d.day, d.year, d.hour, d.minute, d.second))
        w('  Ion Data: Name,   Mass,   Energy (keV)')
        w('           ' + symbol2(ti.ion_z) + '    ' + vb_format(t.mass[0], '###.#####', single=True) +
          ' ' + vb_format(ti.energy, '#########.##'))
        if ti.cascade != 3:
            for i in range(1, t.nel + 1):
                sym = symbol2(t.z[i])
                w('  Displacement Energy of ' + sym + ' = ' + vb_format(t.ed[i], '###0.00000') + ' eV')
                w('  Latt.Binding Energy of ' + sym + ' = ' + vb_format(t.lattice[i], '###0.00000') + ' eV')
                w('  SurfaceBind. Energy of ' + sym + ' = ' + vb_format(t.surface[i], '###0.00000') + ' eV')
        w(self.header_lines[0])
        for l in self.header_lines[2:]:
            w(l)

    def _range3d_header(self):
        ti = self.ti
        t = self.eng.t
        self._open('r3d', 'RANGE_3D.txt')
        w = lambda s='': self._w('r3d', s)
        w(self.g_b2c)
        w(self.g_b30)
        w('====== SRIM Outputs\\RANGE_3D.txt : File of Final Ion Positions =====')
        w('')
        line = ('Ion = ' + symbol2(ti.ion_z) + ' (' + vb_format(ti.ion_z, '00') + ')   ' +
                'Ion Mass= ' + vb_format(t.mass[0], '000.0000', single=True))
        w(line)
        w('Energy  = ' + vb_format(ti.energy, '0000000.E+00') + ' keV')
        w('Ion Angle to Surface = ' + vb_format(ti.angle, '00.0') + ' degrees')
        w('==================== TARGET =======================')
        for j, l in enumerate(ti.layers):
            n = vb_str(j + 1) + ' '
            w('  Layer #' + n + '- ' + l.name)
            w('  Layer #' + n + '- Bottom Depth=' + vb_format(t.depth[j], '000000.E+00') + ' A')
            w('  Layer #' + n + '- Density =' + vb_format(l.density, '0000.0000', single=True) + ' g/cm3')
            for i in range(1, t.nel + 1):
                pct = int(round(f32(t.frac[j][i]) * 100))
                w('  Layer #' + n + '- ' + vb_str(pct) + ' ' + '% of ' + symbol2(t.z[i]))
            w('\xcd' * 52)
        w('                                            ')
        w('Ion       Depth  X   Lateral Y   Lateral Z  ')
        w('Number    (Angstrom) (Angstrom)  (Angstrom) ')
        w('-------  ----------- ----------- -----------')

    # ------------------------------------------------------- per-ion lines
    @staticmethod
    def atom_line(tag, ion, z, e, x, y, zz, cx, cy, cz, single=False):
        """0x43ddd0: one line of TRANSMIT/BACKSCAT/SPUTTER/TRIMOUT"""
        head = (tag + ('     ' + vb_str(ion))[-5:])[:6]
        s = head + ' ' + (' ' + str(z))[-2:] + ' '
        s += vb_format(e, '.0000000E+00 ', single=single)
        s += vb_format(x, '  ######0E+00 ', single=single)
        s += vb_format(y, ' .0000E+00 ;-.0000E+00 ', single=single)
        s += vb_format(zz, ' .0000E+00;-.0000E+00', single=single)
        s += vb_format(cx, '   .0000000 ;  -.0000000 ')
        s += vb_format(cy, ' .0000000 ;-.0000000 ')
        s += vb_format(cz, ' .0000000;-.0000000')
        return s

    def transmitted(self, ion, st):
        """0x440300"""
        if 'trans' not in self.files and 'out' not in self.files:
            return
        line = self.atom_line('T', ion, self.ti.ion_z, st.e, st.x, st.y, st.z,
                              cos_string(st.cx, -0.999999), cos_string(st.cy), cos_string(st.cz))
        if 'out' in self.files:
            self._w('out', line)
        if 'trans' in self.files:
            self._w('trans', line)

    def backscattered(self, ion, st):
        """0x43d850"""
        line = self.atom_line('B', ion, self.ti.ion_z, st.e, st.x, st.y, st.z,
                              cos_string(st.cx, -0.999999), cos_string(st.cy), cos_string(st.cz))
        if 'back' in self.files:
            self._w('back', line)
        if 'out' in self.files:
            self._w('out', line)

    def sputtered(self, ion, z, e, x, y, zz, cx, cy, cz):
        """0x441ff0 (all arguments are TRIM Singles)"""
        def clamp(c):
            if c > 0.9999999:
                return 0.9999999
            if c < -0.9999999:
                return -0.9999999
            return c
        if x < 0:
            x = 0
        line = self.atom_line('S', ion, z, e, x, y, zz, clamp(cx), clamp(cy), clamp(cz), single=True)
        if 'sput' in self.files:
            self._w('sput', line)
        if 'out' in self.files:
            self._w('out', line)

    def range3d(self, ion, st):
        """0x43d310"""
        if 'r3d' not in self.files:
            return
        self._w('r3d', vb_format(ion, '0000000') + vb_format(st.x, '   0.0000E+00') +
                vb_format(st.y, '  0.0000E+00; -0.0000E+00') + vb_format(st.z, '  0.0000E+00; -0.0000E+00'))


# ===================================================================
# distribution files written at the end of the calculation (0x46a190)
# ===================================================================
_sqrt = math.sqrt
DASH68 = '-' * 68
EQ68 = '=' * 68


class TrimDistributions(object):
    """RANGE.txt, VACANCY.txt, IONIZ.txt, PHONON.txt, E2RECOIL.txt,
    NOVAC.txt and LATERAL.txt (TRIM's "save all plot data" routine)."""

    def __init__(self, engine, directory):
        self.eng = engine
        self.ti = engine.ti
        self.t = engine.t
        self.tal = engine.tal
        self.dir = directory
        ti, t, tal = self.ti, self.t, self.tal
        eng = engine
        nel, nlay = t.nel, t.nlay
        self.binw = eng.binw
        # depth of the bins (+0.01 A)
        self.depth = [None] + [(eng.xmin + eng.binw * i) + 0.01 for i in range(1, 101)]
        # number of ions (minus the remaining fraction of the last one)
        e0kev = ti.energy
        self.ions = eng.ion_number - (eng.last_e / 1000) / e0kev
        if self.ions <= 0:
            self.ions = 1
        self.kp = eng.kpflag
        if self.kp == 1:
            self.kp_lines = [DASH68,
                             '    Recoil/Damage Calculations made with Kinchin-Pease Estimates    ',
                             DASH68]
        else:
            self.kp_lines = []
        # target description lines (0x4fe5ac per layer, 0x4fe5a0 per element)
        self.meanmass = []
        for j in range(nlay):
            m = 0
            for i in range(1, nel + 1):
                m = m + t.frac[j][i] * t.mass[i]
            self.meanmass.append(m)
        self.dens_line = []
        self.elem_lines = []
        for j, l in enumerate(ti.layers):
            self.dens_line.append(self._dens_line(j))
            el = []
            for i in range(1, nel + 1):
                fr = t.frac[j][i]
                if fr == 0:
                    el.append('')
                    continue
                s = '  Layer #' + vb_str(j + 1) + '- ' + symbol2(t.z[i]) + ' ='
                s += (vb_str(f32(fr * 100), True) + '    ')[:5] + ' Atomic Percent ='
                mf = (fr * t.mass[i]) / self.meanmass[j]
                s += (vb_str(mf * 100) + '    ')[:5] + ' Mass Percent'
                el.append(s)
            self.elem_lines.append(el)
        self.files = {}

    # -------------------------------------------------------------- helpers
    def _dens_line(self, j, nchar=6):
        """'  Layer # 1- Density = 4.976E22 atoms/cm3 = 2.321 g/cm3' (TDATA
        keeps 11 characters of the mass density when it is below 0.001)"""
        n = self.t.N[j]
        if n >= 0.01:
            dl = vb_str(n * 100)[:6] + 'E22 atoms/cm3'
        else:
            dl = vb_str(n / 10)[:6] + 'E19 atoms/cm3'
        line = '  Layer #' + vb_str(j + 1) + '- Density =' + dl
        return line + ' =' + vb_str(self.ti.layers[j].density, True)[:nchar] + ' g/cm3'

    def _open(self, name):
        f = open(os.path.join(self.dir, name), 'w', encoding='latin-1', newline='')
        lines = []
        return f, lines

    @staticmethod
    def _save(f, lines):
        f.write(''.join(l + CRLF for l in lines))
        f.close()

    def _title_block(self, out, title2):
        ti = self.ti
        out.append('====== ' + ti.target_name + ' =======')
        out.append('        ' + 'SRIM-2013.00')
        for l in title2:
            out.append(l)
        out.append(self.eng.out.title if self.eng.out is not None else self._title())
        out.append('-' * 78)
        out.extend(self.kp_lines)

    def _title(self):
        ti = self.ti
        names = '+'.join(l.name for l in ti.layers)
        title = ('====== TRIM Calc.= ' + symbol2(ti.ion_z) + '(' + energy_label(ti.energy * 1000) +
                 ') ==> ' + names[:24] + '(' + depth_label(self.t.depth[-1]) + ') ')
        return (title + '=' * 78)[:78]

    def _ion_line(self):
        # RANGE.txt writes the energy with its leading zero (0.1), the other
        # files use Str$(E) (.1)
        e = vb_cstr(self.ti.energy)
        if e.startswith('.'):
            e = '0' + e
        return ' Ion    = ' + symbol2(self.ti.ion_z) + '   Energy = ' + e + ' ' + ' keV'

    def _target_block(self, out, width_sep):
        """layers, widths, densities and composition"""
        ti, t = self.ti, self.t
        out.append('============= TARGET MATERIAL ======================================')
        for j, l in enumerate(ti.layers):
            out.append('Layer ' + ('  ' + str(j + 1))[-2:] + ' : ' + l.name)
            w = vb_format(t.width[j], '0.E+00')
            w = ('%10s' % w)[-10:]
            if width_sep:
                out.append('Layer Width = ' + w + ' A ' + self.dens_line[j])
            else:
                out.append('Layer Width = ' + w + ' A ;')
                out.append(self.dens_line[j])
            for s in self.elem_lines[j]:
                if s:
                    out.append(s)
        out.append(EQ68)

    # ---------------------------------------------------------------- RANGE
    def write_range(self):
        """0x46ebe0: RANGE.txt"""
        ti, t, tal = self.ti, self.t, self.tal
        f, out = self._open('RANGE.txt')
        self._title_block(out, ['==================================================',
                                '      ION and final RECOIL ATOM Distributions     ',
                                ' See  SRIM Outputs\\TDATA.txt for calculation details',
                                '=================================================='])
        if self.kp == 1:
            out.append('-' * 64)
            out.append('  There are NO RECOILS unless a full TRIM calculation is made.   ')
            out.append('-' * 64)
        out.append('See file :  SRIM Outputs\\TDATA.txt   for details of calculation')
        out.append(self._ion_line())
        self._target_block(out, True)
        out.append(' Total Ions calculated =' + vb_format(self.ions, '######.00'))
        f11 = lambda v: ('%11s' % vb_format(v, '000.0E+00'))[-11:]
        out.append(' Ion Average Range = ' + f11(tal.g_a38) + ' A   ' + 'Straggling = ' + f11(tal.g_a48) + ' A')
        out.append(' Ion Lateral Range = ' + f11(tal.g_a10) + ' A   ' + 'Straggling = ' + f11(tal.g_a20) + ' A')
        out.append(' Ion Radial  Range = ' + f11(tal.g_9f8) + ' A   ' + 'Straggling = ' + f11(tal.g_a08) + ' A')
        n = (self.ions - tal.n_back) - tal.n_trans
        if n > 3:
            out.append(EQ68)
            if tal.n_trans + tal.n_back > 0:
                out.append(' Transmitted Ions =' + vb_format(tal.n_trans, '######') +
                           ';  Backscattered Ions =' + vb_format(tal.n_back, '#####'))
                out.append(' (These are not included in Skewne- and Kurtosis below.)')
                out.append('')
            avg = tal.g280 / n
            tal.g_a38 = avg
            strag = _sqrt(abs(tal.g288 / n - avg ** 2))
            if not (strag > 0):
                strag = 1e-7
            tal.g_a48 = strag
            s1, s2, s3, s4 = tal.g280, tal.g288, tal.g290, tal.g298
            a = avg
            x = ((s3 - (a * 3.0) * s2) + ((a ** 2.0) * 3.0) * s1)
            x = x - (a ** 3.0) * n
            skew = x / ((strag ** 3.0) * n)
            out.append(' Range Skewne- = ' + vb_format(skew, '000.0000') + '  &' + '= [-(X-Rp)^3]/[N*Straggle^3]')
            y = ((s4 - (a * 4.0) * s3) + ((a ** 2.0) * 6.0) * s2)
            y = y - ((a ** 3.0) * 4.0) * s1
            y = y + (a ** 4.0) * n
            kurt = y / ((strag ** 4.0) * n)
            out.append(' Range Kurtosis = ' + vb_format(kurt, '000.0000') + '  &' + '= [-(X-Rp)^4]/[N*Straggle^4]')
            out.append(' Statistical definitions above are those used in VLSI implant modelling.')
            out.append(EQ68)
        out.append('=================================================================')
        out.append('  Table Distribution Units are >>>  (Atoms/cm3) / (Atoms/cm2)  <<<  ')
        out.append('=================================================================')
        h1 = '   DEPTH         ' + symbol2(ti.ion_z) + '    '
        h2 = '   (Ang.)       Ions   '
        h3 = '-----------  ----------'
        if self.kp == 0:
            for i in range(1, t.nel + 1):
                h1 = h1 + '     ' + symbol2(t.z[i]) + '    '
                h2 = h2 + '  Tgt. Atoms'
                h3 = h3 + '  ----------'
        if self.kp > 0:
            h1 = h1 + '   Recoil     '
            h2 = h2 + '  Distribution'
            h3 = h3 + '  ------------'
        out.extend([h1, h2, h3])
        ncol = t.nel + 1 if self.kp == 0 else 2
        for i in range(1, 101):
            line = vb_format(self.depth[i], '000000.E+00')
            for j in range(1, ncol + 1):
                v = (1e8 * tal.d98c.get(j, i)) / self.ions
                line += vb_format(v, '  0.0000E+00')
            out.append(line)
        self._save(f, out)

    # ------------------------------------------------------ shared pieces
    def _head(self, out, lines, see):
        self._title_block(out, lines)
        if see is None:
            return
        out.append(see)
        out.append(' Ion    = ' + symbol2(self.ti.ion_z) + '   Energy =' + vb_str(self.ti.energy) + ' keV')

    def _target_block2(self, out, pad_layer, width_fmt, width_rset, unit=' A ;'):
        """target description used by IONIZ/VACANCY/PHONON/E2RECOIL/NOVAC"""
        ti, t = self.ti, self.t
        out.append('============= TARGET MATERIAL ======================================')
        for j, l in enumerate(ti.layers):
            if pad_layer:
                out.append('Layer ' + ('  ' + str(j + 1))[-2:] + ' : ' + l.name)
            else:
                out.append('Layer ' + str(j + 1) + ' : ' + l.name)
            w = vb_format(t.width[j], width_fmt)
            if width_rset:
                w = ('%10s' % w)[-10:]
            out.append('Layer Width = ' + w + unit)
            out.append(self.dens_line[j])
            for s in self.elem_lines[j]:
                if s:
                    out.append(s)
        out.append(EQ68)

    def _totals(self, out, note=()):
        tal = self.tal
        inv = self.eng.invbinw
        vac = int(round((tal.g244 / inv) / self.ions))
        out.append(' Total Target Vacancies     =' + vb_str(vac) + ' ' + '/Ion')
        out.extend(note)
        if self.kp != 1:
            disp = int(round(((tal.g244 + tal.g1d8) / inv) / self.ions))
            out.append(' Total Target Displacements =' + vb_str(disp) + ' ' + '/Ion')
            rep = int(round((tal.g1d8 / inv) / self.ions))
            out.append(' Total Target Replacement Collisions =' + vb_str(rep) + ' ' + '/Ion')

    # ---------------------------------------------------------------- IONIZ
    def write_ioniz(self):
        """0x470a20: IONIZ.txt"""
        tal = self.tal
        f, out = self._open('IONIZ.txt')
        self._head(out, ['========================================',
                         '      Ion and Recoil IONIZATION     ',
                         ' See SRIM Outputs\\TDATA.txt for details   ',
                         '========================================'],
                   ' See file :  SRIM Outputs\\TDATA.txt   for calculation data')
        self._target_block2(out, False, '00000.E+00', False)
        out.append(' Total Ions calculated =' + vb_format(self.ions, '000000.00'))
        out.append('=' * 63)
        out.append('  Ionization Energy Units are  >>>>  eV /(Angstrom-Ion)  <<<<  ')
        out.append('=' * 63)
        if self.ions < 20:
            out.append('-' * 63)
            out.append('-  NOTE: Because of the application of straggling to the ion  -')
            out.append('-  energy, Ionization is not accurate for less than 20 ions!  -')
            out.append('-' * 63)
        out.append('')
        out.append('  TARGET       IONIZ.       IONIZ.  ')
        out.append('  DEPTH         by           by     ')
        out.append('  (Ang)        IONS         RECOILS ')
        out.append('-----------  -----------  -----------')
        for i in range(1, 101):
            out.append(vb_format(self.depth[i], '000000.E+00') +
                       vb_format(tal.d964.get(1, i) / self.ions, '  0000.00E+00') +
                       vb_format(tal.d964.get(2, i) / self.ions, '  0000.00E+00'))
        self._save(f, out)

    # -------------------------------------------------------------- VACANCY
    def write_vacancy(self):
        """0x471770: VACANCY.txt"""
        t, tal = self.t, self.tal
        f, out = self._open('VACANCY.txt')
        self._head(out, ['==========================================',
                         '    Ion and Target VACANCY production     ',
                         '  See SRIM Outputs\\TDATA.txt for calc. details',
                         '=========================================='],
                   ' See file :  SRIM Outputs\\TDATA.txt for calculation data')
        self._target_block2(out, True, '0.E+00', True)
        out.append(' Total Ions calculated =' + vb_format(self.ions, '0.00'))
        note = []
        if self.kp == 2:
            note = ['=' * 73,
                    '  NOTE: Vacancies not accurately calculated in SPUTTERING calculation.   ',
                    '  Atoms move if Energy > E(surface) instead of Energy > E(displacement)  ',
                    '  Quick Fix for near-surface vacancies: Make E(surface)= E(displacement),',
                    '      and use Monolayer Steps (TRIM Damage Setup). Not perfect, but OK.  ',
                    '=' * 73]
        self._totals(out, note)
        out.append('')
        if self.kp == 0:
            out.append('!!!! NOTE : 2nd Column below is number of Primary Knock-Ons !!!!')
            out.append('           ( PKO are number of Target Atoms Recoiling from the Ion. )')
        out.append('=' * 58)
        out.append('  Table Units are  >>>> Vacancies/(Angstrom-Ion)  <<<<  ')
        out.append('=' * 58)
        h = ['   TARGET    ', '   DEPTH     ', '   (Ang.)    ', '-----------  ']
        if self.kp == 0:
            for i in range(0, t.nel + 1):
                h[1] = h[1] + '    ' + symbol2(t.z[i]) + '       '
            h[2] = h[2] + ' Knock-Ons   '
            for i in range(1, t.nel + 1):
                h[2] = h[2] + ' Vacancies   '
            for i in range(0, t.nel + 1):
                h[3] = h[3] + '-----------  '
        fmt = '0000.00E+00  '
        if self.ti.cascade == 3 or self.kp == 1:
            h[0] = h[0] + ' VACANCIES     VACANCIES '
            h[1] = h[1] + '    by            by     '
            h[2] = h[2] + '   IONS         RECOILS  '
            h[3] = h[3] + '-----------  ------------'
        out.extend(h)
        ncol = t.nel + 1 if self.kp == 0 else 2
        if self.ti.cascade == 3:
            ncol = 2
        for i in range(1, 101):
            line = ('%13s' % vb_format(self.depth[i], '000000.E+00  '))
            for j in range(1, ncol + 1):
                line += vb_format(tal.d9a0.get(j, i) / self.ions, fmt)
            out.append(line)
        out.append('')
        out.append(' To convert to Energy Lost - multiply by Average Binding Energy = ' +
                   vb_str(self.eng.avg_lattice, isinstance(self.eng.avg_lattice, F32)) + ' ' + ' eV/Vacancy')
        self._save(f, out)

    # --------------------------------------------------------------- PHONON
    def write_phonon(self):
        """0x472de0: PHONON.txt"""
        tal = self.tal
        f, out = self._open('PHONON.txt')
        self._head(out, ['======================================== ',
                         '       Ion and Recoil PHONONS        ',
                         ' See SRIM Outputs\\TDATA.txt for calc. details',
                         '======================================== '],
                   ' See file :  SRIM Outputs\\TDATA.txt calculation data')
        self._target_block2(out, True, '0.E+00', True, ' AAA ;')
        out.append(' Total Ions calculated =' + vb_format(self.ions, '0.00'))
        out.append('=' * 56)
        out.append('  Table Units are  >>>> Phonons /(Angstrom-Ion)  <<<<  ')
        out.append('=' * 56)
        out.append('  DEPTH       PHONONS     PHONONS    ')
        out.append('  (Ang.)      by ION      by RECOILS ')
        out.append('-----------  -----------  -----------')
        r = lambda v: ('%11s' % vb_format(v, '000000.E+00'))[-11:]
        for i in range(1, 101):
            out.append(r(self.depth[i]) + '  ' + r(tal.d970.get(1, i) / self.ions) + '  ' +
                       r(tal.d970.get(2, i) / self.ions))
        self._save(f, out)

    # ------------------------------------------------------------- E2RECOIL
    def write_e2recoil(self):
        """0x473b00: E2RECOIL.txt"""
        t, tal = self.t, self.tal
        f, out = self._open('E2RECOIL.txt')
        self._head(out, ['=======================================',
                         '    Energy transferred to RECOILS      ',
                         ' See SRIM Outputs\\TDATA.txt for calc. details ',
                         '======================================='],
                   ' See file: SRIM Outputs\\TDATA.txt for data about the calculation')
        self._target_block2(out, True, '00000.E+00', False, ' AAA ;')
        out.append(' Total Ions calculated =' + vb_format(self.ions, '0.00'))
        out.append('=' * 59)
        out.append('  Table Energy Units are  >>>>  eV /(Angstrom-Ion)  <<<<  ')
        out.append('=' * 59)
        h = ['   DEPTH       ENERGY      ', '   DEPTH        from       ', '   (Ang.)       IONS       ',
             '----------   ------------  ']
        fmt = '0000.00E+00  '
        if self.kp == 0:
            for i in range(1, t.nel + 1):
                h[0] += 'Energy       '
                h[1] += 'Absorbed     '
                h[2] += 'by  ' + symbol2(t.z[i]) + '       '
                h[3] += '-----------  '
        w2 = 13
        if self.kp == 2:
            # TRIM leaves the column format and field of the previously
            # written files (RANGE: '  0.0000E+00', LATERAL: 14 characters)
            fmt = '  0.0000E+00'
            w2 = 14
        if self.kp == 1:
            h[0] += ' ENERGY    '
            h[1] += ' Absorbed  '
            h[2] += ' by Recoils'
            h[3] += '-----------'
            fmt = '0000.00E+00'
            w2 = 11
        out.extend(h)
        ncol = t.nel + 1 if self.kp == 0 else 2
        for i in range(1, 101):
            line = ('%13s' % vb_format(self.depth[i], '000000.E+00  '))[-13:]
            # RSet: a negative value (round-off of the energy balance) is
            # one character too long and loses its rightmost character
            line += _rset(vb_format(tal.d95c.get(1, i) / self.ions, '0000.00E+00  '), 13)
            for j in range(2, ncol + 1):
                line += _rset(vb_format(tal.d95c.get(j, i) / self.ions, fmt), w2)
            out.append(line)
        self._save(f, out)

    # ---------------------------------------------------------------- NOVAC
    def write_novac(self):
        """0x474c30: NOVAC.txt"""
        tal = self.tal
        f, out = self._open('NOVAC.txt')
        self._head(out, ['=========================================',
                         '        Replacement Collisions       ',
                         ' See SRIM Outputs\\TDATA.txt for calc. details   ',
                         '========================================='],
                   None)
        if self.kp == 1:
            out.append('')
            out.append(' There are NO Replacement Collisions calculated.')
            out.append(' You must execute a full TRIM cascade calculation for these.')
            self._save(f, out)
            return
        out.append(' See file: SRIM Outputs\\TDATA.txt for data about the calculation')
        out.append(' Ion    = ' + symbol2(self.ti.ion_z) + '   Energy =' + vb_str(self.ti.energy) + ' keV')
        self._target_block2(out, True, '0.E+00', True, ' AAA ;')
        out.append(' Total Ions calculated =' + vb_format(self.ions, '0.00'))
        self._totals(out)
        out.append('=' * 64)
        out.append('   Table Vacancy Units are  >>>>  Number /(Angstrom-Ion)  <<<<  ')
        out.append('=' * 64)
        out.append('  DEPTH (A)      Number')
        out.append('------------   -----------')
        for i in range(1, 101):
            out.append(('%15s' % vb_format(self.depth[i], '0000000.E+00   '))[-15:] +
                       ('%11s' % vb_format(tal.d96c.get(i) / self.ions, '0000.000000'))[-11:])
        self._save(f, out)

    # -------------------------------------------------------------- LATERAL
    def write_lateral(self):
        """0x4759f0: LATERAL.txt"""
        ti, t, tal = self.ti, self.t, self.tal
        f, out = self._open('LATERAL.txt')
        self._title_block(out, ['=====================================',
                                '       LATERAL SPREAD OF IONS        ',
                                ' See SRIM Outputs\\TDATA.txt for calc. details   ',
                                '\xcc===================================\xca' + '=' * 35,
                                '  NOTE : The Projected LATERAL Spread is based on a projection on a     ',
                                '  plane- just the way the XY or XZ displays look on the screen.  The    ',
                                '  Lateral Projected Range is defined as the average of the absolute     ',
                                '  values of the projected lateral displacements from the X-axis.        ',
                                '  One use of this range: If a mask covers half of a target, the lateral ',
                                '  range is the distance under the mask where the concentration drops to ',
                                '  25% (at the mask edge the concentration= 50% of the uniform value).   ',
                                ' ' * 72,
                                '  The Radial Range is the mean radial displacement range from the       ',
                                '  X-axis a-uming cylindrical symmetry.                                 ',
                                '=' * 72])
        out.append(' See file: SRIM Outputs\\TDATA.txt for data about the calculation')
        out.append(' Ion    = ' + symbol2(ti.ion_z) + '   Energy =' + vb_str(ti.energy) + ' keV')
        out.append(' Total Ions calculated =' + vb_format(self.ions, '0.00'))
        self._target_block2(out, True, '00000.E+00', False, ' AAA ;')
        f11 = lambda v: ('%11s' % vb_format(v, '000000.E+00'))[-11:]
        out.append(' Ion Average Range = ' + f11(tal.g_a38) + ' A   ' + 'Straggling = ' + f11(tal.g_a48) + ' A')
        out.append(' Ion Lateral Range = ' + f11(tal.g_a10) + ' A   ' + 'Straggling = ' + f11(tal.g_a20) + ' A')
        out.append(' Ion Radial  Range = ' + f11(tal.g_9f8) + ' A   ' + 'Straggling = ' + f11(tal.g_a08) + ' A')
        out.append(EQ68)
        out.append('                LATERAL     Projected      LATERAL        Radial    ')
        out.append('  TARGET       PROJ.RANGE   Straggling      RADIAL      Straggling  ')
        out.append('  DEPTH (A)   (Angstroms)   (Angstroms)   (Angstroms)   (Angstroms) ')
        out.append('------------ ------------- ------------- ------------- -------------')
        r14 = lambda v: ('%14s' % vb_format(v, '00000000.E+00 '))[-14:]
        inv = self.eng.invbinw
        for i in range(1, 101):
            n = tal.d98c.get(1, i) / inv
            if n == 0:
                a = c = b = d = 0
            else:
                a = tal.d9b0.get(1, i) / n
                b = tal.d9ac.get(1, i) / n
                c = _sqrt(tal.d9b0.get(2, i) / n)
                d = _sqrt(abs(tal.d9ac.get(2, i) / n - (tal.d9ac.get(1, i) / n) ** 2))
            out.append(('%13s' % vb_format(self.depth[i], '0000000.E+00 '))[-13:] +
                       r14(a) + r14(c) + r14(b) + r14(d))
        self._save(f, out)

    # ---------------------------------------------------------------- TDATA
    def write_tdata(self):
        """0x46c090: TDATA.txt"""
        ti, t, tal, eng = self.ti, self.t, self.tal, self.eng
        f, out = self._open('TDATA.txt')
        self._title_block(out, ['=======================================',
                                '      Data  of  TRIM  Calculation      ',
                                '======================================='])
        out.append('Ion = ' + symbol2(ti.ion_z) + ' (' + (' ' + str(ti.ion_z))[-2:] + ')   ' +
                   'Ion Mass= ' + vb_format(ti.ion_mass, '000.0000', single=True))
        out.append('Energy  = ' + ('%12s' % vb_format(ti.energy, '0000000.E+00'))[-12:] + ' keV')
        out.append('Ion Angle to Surface = ' + vb_cstr(ti.angle) + ' degrees')
        out.append('============= TARGET MATERIAL ======================================')
        r11 = lambda v, fmt: ('%11s' % vb_format(v, fmt))[-11:]
        for j, l in enumerate(ti.layers):
            out.append('  Layer #' + vb_print_num(j + 1) + '- ' + l.name)
            out.append('  Layer #' + vb_print_num(j + 1) + '- Bottom Depth=' + r11(t.depth[j], '0.E+00') + ' A')
            out.append(self._dens_line(j, 11 if f32(0.001) > l.density else 6))
            for s in self.elem_lines[j]:
                if s:
                    out.append(s)
            out.append(EQ68)
        for i in range(1, t.nel + 1):
            out.append('Target energies for target atom = ' + symbol2(t.z[i]))
            out.append('    Displacement =' + vb_print_num(t.ed[i]) + 'eV, ' +
                       'Binding =' + vb_print_num(t.lattice[i]) + 'eV, ' +
                       'Surface =' + vb_print_num(t.surface[i]) + 'eV')
        out.append(EQ68)
        out.append('Depth Range of Tabulated Data=  ' + r11(eng.xmin, '000000.E+00') + '  -  ' +
                   r11(eng.xmax, '000000.E+00') + '  Angstroms')
        out.append(EQ68)
        out.append('Total Ions calculated =  ' + vb_format(self.ions, '000000'))
        out.append('Average Range         =  ' + r11(tal.g_a38, '00000.E+00') + ' Angstroms')
        out.append('Average Straggling    =  ' + r11(tal.g_a48, '00000.E+00') + ' Angstroms')
        out.append('Average Vacancy/Ion   =  ' + r11(tal.g_a50, '00000.E+00'))
        out.append(EQ68)
        out.append('Total Backscattered Ions= ' + vb_print_num(tal.n_back))
        out.append('Total Transmitted Ions  = ' + vb_print_num(tal.n_trans))
        out.append(EQ68)
        out.append('See files SRIM Outputs\\*.txt for tables :')
        out.append('   TDATA.txt    -- Details of the calculation.')
        out.append('   RANGE.txt    -- Final distribution of Ions/Recoils')
        out.append('   VACANCY.txt  -- Distribution of Vacancies')
        out.append('   IONIZ.txt    -- Distribution of Ionization')
        out.append('   PHONON.txt   -- Distribution of Phonons')
        out.append('   E2RECOIL.txt -- Energy transferred to recoils')
        out.append('   NOVAC.txt    -- Replacement collisions')
        out.append('   LATERAL.txt  -- Lateral Spread of Ions')
        self._save(f, out)

    def write_all(self):
        """all distribution files, in TRIM's order (0x46be80)"""
        self.write_range()
        self.write_lateral()
        self.write_ioniz()
        self.write_phonon()
        self.write_e2recoil()
        self.write_vacancy()
        self.write_novac()
        self.write_tdata()
