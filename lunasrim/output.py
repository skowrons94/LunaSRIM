"""Output files of TRIM (SRIM-2013 formats).

The text written here reproduces the files TRIM.exe produces in its
``SRIM Outputs`` directory.  Numbers are formatted with the Visual Basic
conventions (see :mod:`lunasrim.vbfmt`).
"""
import datetime
import os

from . import data_files
from .trimin import f32
from .vbfmt import vb_format, vb_str, vb_cstr

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
