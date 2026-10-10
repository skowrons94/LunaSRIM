"""Readers for the TRIM output files (``SRIM Outputs`` directory).

They work for files written by lunasrim as well as by SRIM-2013 itself and
follow the attribute names of pysrim's ``srim.output`` module.  Columns are
returned as numpy arrays when numpy is installed, as lists otherwise.
"""
import os
import re

try:
    import numpy as _np
except ImportError:          # pragma: no cover
    _np = None


def _array(v):
    return _np.array(v, dtype=float) if _np is not None else list(v)


def _read(path):
    with open(path, 'rb') as fh:
        return fh.read().decode('latin-1').replace('\r\n', '\n').split('\n')


def _floats(line):
    try:
        return [float(x) for x in line.split()]
    except ValueError:
        return None


def _table(lines):
    """rows of the numeric table that follows the last '-----' ruler"""
    start = None
    for i, l in enumerate(lines):
        s = l.strip()
        if s and set(s) <= set('- ') and s.startswith('---'):
            start = i + 1
    rows = []
    if start is None:
        return rows
    for l in lines[start:]:
        v = _floats(l)
        if not v:
            if rows:
                break
            continue
        rows.append(v)
    return rows


def _columns(rows, n=None):
    if not rows:
        return []
    n = n or max(len(r) for r in rows)
    return [_array([r[k] if k < len(r) else 0.0 for r in rows]) for k in range(n)]


def _header_value(lines, label, cast=float):
    for l in lines:
        if label in l:
            m = re.search(re.escape(label) + r'\s*=?\s*([-+0-9.Ee]+)', l)
            if m:
                try:
                    return cast(m.group(1))
                except ValueError:
                    return None
    return None


class _TableFile(object):
    name = None

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        self.lines = _read(self.path)
        self.data = _columns(_table(self.lines))
        self.num_ions = _header_value(self.lines, 'Total Ions calculated')
        self.depth = self.data[0] if self.data else _array([])


class Ioniz(_TableFile):
    """IONIZ.txt: ionization energy loss (eV/(Angstrom-Ion)) by ions and recoils"""
    name = 'IONIZ.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.ions = self.data[1]
        self.recoils = self.data[2]


class Phonons(_TableFile):
    """PHONON.txt: phonons (eV/(Angstrom-Ion)) produced by ions and recoils"""
    name = 'PHONON.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.ions = self.data[1]
        self.recoils = self.data[2]


class EnergyToRecoils(_TableFile):
    """E2RECOIL.txt: energy transferred to recoils (eV/(Angstrom-Ion));
    ``absorbed`` has one column per target element (full cascades) or the
    total (Kinchin-Pease)"""
    name = 'E2RECOIL.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.ions = self.data[1]
        self.absorbed = self.data[2:]


class Vacancy(_TableFile):
    """VACANCY.txt: vacancies/(Angstrom-Ion).  Full cascade: ``knock_ons``
    plus one ``vacancies`` column per element; Kinchin-Pease: vacancies by
    ions and by recoils"""
    name = 'VACANCY.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.knock_ons = self.data[1]
        self.vacancies = self.data[2:]
        self.vacancies_per_ion = _header_value(self.lines, 'Total Target Vacancies')
        self.displacements_per_ion = _header_value(self.lines, 'Total Target Displacements')
        self.replacements_per_ion = _header_value(self.lines, 'Total Target Replacement Collisions')


class NoVacancy(_TableFile):
    """NOVAC.txt: replacement collisions/(Angstrom-Ion) (full cascades only)"""
    name = 'NOVAC.txt'

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        self.lines = _read(self.path)
        self.num_ions = _header_value(self.lines, 'Total Ions calculated')
        rows = [] if any('NO Replacement' in l for l in self.lines) else _table(self.lines)
        self.data = _columns(rows, 2)
        self.depth = self.data[0] if rows else _array([])
        self.number = self.data[1] if rows else _array([])


class Range(_TableFile):
    """RANGE.txt: final ion (and recoil) distributions,
    (atoms/cm3)/(atoms/cm2)"""
    name = 'RANGE.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.ions = self.data[1]
        self.elements = self.data[2:]
        self.average_range = _header_value(self.lines, 'Ion Average Range')
        self.lateral_range = _header_value(self.lines, 'Ion Lateral Range')
        self.radial_range = _header_value(self.lines, 'Ion Radial  Range')
        m = [re.findall(r'Straggling =\s*([-+0-9.Ee]+)', l) for l in self.lines]
        s = [float(x[0]) for x in m if x]
        self.straggling = s[0] if len(s) > 0 else None
        self.lateral_straggling = s[1] if len(s) > 1 else None
        self.radial_straggling = s[2] if len(s) > 2 else None


class Lateral(_TableFile):
    """LATERAL.txt: lateral spread of the ions versus depth (Angstrom)"""
    name = 'LATERAL.txt'

    def __init__(self, directory, filename=None):
        _TableFile.__init__(self, directory, filename)
        self.lateral_range = self.data[1]
        self.lateral_straggling = self.data[2]
        self.radial_range = self.data[3]
        self.radial_straggling = self.data[4]


class _AtomFile(object):
    """TRANSMIT.txt / BACKSCAT.txt / SPUTTER.txt / TRIMOUT.txt"""
    name = None
    tags = 'TBS'

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        rows = []
        kinds = []
        for l in _read(self.path):
            if len(l) > 2 and l[0] in self.tags and l[1] == ' ':
                v = _floats(l[1:])
                if v and len(v) == 9:
                    rows.append(v)
                    kinds.append(l[0])
        self.kind = kinds
        cols = _columns(rows, 9) if rows else [_array([]) for k in range(9)]
        self.ion_number = cols[0]
        self.z = cols[1]
        self.energy = cols[2]            # eV
        self.position = [cols[3], cols[4], cols[5]]   # x (depth), y, z in Angstrom
        self.direction = [cols[6], cols[7], cols[8]]  # cosines

    def __len__(self):
        return len(self.kind)


class Transmit(_AtomFile):
    name = 'TRANSMIT.txt'
    tags = 'T'


class Backscat(_AtomFile):
    name = 'BACKSCAT.txt'
    tags = 'B'


class Sputter(_AtomFile):
    name = 'SPUTTER.txt'
    tags = 'S'


class TrimOut(_AtomFile):
    name = 'TRIMOUT.txt'


class Range3D(object):
    """RANGE_3D.txt: final position of every stopped ion"""
    name = 'RANGE_3D.txt'

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        rows = [r for r in _table(_read(self.path)) if len(r) == 4]
        cols = _columns(rows, 4) if rows else [_array([]) for k in range(4)]
        self.ion_number = cols[0]
        self.positions = [cols[1], cols[2], cols[3]]

    def __len__(self):
        return len(self.ion_number)


def _num(field):
    """number of a fixed width field; TRIM may put the minus sign of a
    negative value in front of the field's blanks ('-  3.1590E-01')"""
    return float(field.replace(' ', ''))


_EXYZ_NUM = re.compile(r'(-?)\s*(\d\.\d+E[+-]\d+)')


class EXYZ(object):
    """EXYZ.txt: ion energy (keV), position (A), electronic stopping (eV/A)
    and energy of the last recoil (eV) each time the ion has lost the energy
    interval chosen in TRIM.IN"""
    name = 'EXYZ.txt'

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        rows = []
        for l in _read(self.path):
            if l[:7].isdigit():
                v = [float(sg + x) for sg, x in _EXYZ_NUM.findall(l[7:])]
                if len(v) == 6:
                    rows.append([int(l[:7])] + v)
        cols = list(zip(*rows)) if rows else [()] * 7
        self.ion_number = _array(cols[0])
        self.energy = _array(cols[1])
        self.positions = [_array(cols[2]), _array(cols[3]), _array(cols[4])]
        self.se = _array(cols[5])
        self.recoil_energy = _array(cols[6])

    def __len__(self):
        return len(self.ion_number)

    def ion(self, n):
        """rows of ion n as (energy, x, y, z, se, recoil_energy) lists"""
        idx = [i for i, k in enumerate(self.ion_number) if k == n]
        return [[float(c[i]) for i in idx] for c in
                (self.energy, self.positions[0], self.positions[1], self.positions[2], self.se,
                 self.recoil_energy)]


class Collisions(object):
    """COLLISON.txt: ion collisions which displaced a target atom.  For each
    collision: ion number, energy (keV), position (A), electronic stopping
    (eV/A), atom hit, recoil energy (eV) and, depending on the calculation,
    the damage it caused (displacements, vacancies, replacements,
    interstitials; Kinchin-Pease: vacancies only)."""
    name = 'COLLISON.txt'

    def __init__(self, directory, filename=None):
        self.path = os.path.join(directory, filename or self.name)
        self.rows = []
        for l in _read(self.path):
            if not l.startswith('\xb3') or 'Summary' in l:
                continue
            f = l.split('\xb3')[1:-1]
            if len(f) < 8 or not f[0].strip().isdigit():
                continue
            row = {'ion': int(f[0]), 'energy': _num(f[1]), 'x': _num(f[2]), 'y': _num(f[3]),
                   'z': _num(f[4]), 'se': _num(f[5]), 'atom': f[6].strip(), 'recoil_energy': _num(f[7])}
            rest = [x.strip() for x in f[8:] if x.strip()]
            try:
                vals = [float(x) for x in rest]
            except ValueError:
                vals = []           # '<== Start of New Cascade'
            if len(vals) == 4:
                row.update(displacements=vals[0], vacancies=vals[1], replacements=vals[2],
                           interstitials=vals[3])
            elif len(vals) == 1:
                row.update(vacancies=vals[0])
            self.rows.append(row)

    def __len__(self):
        return len(self.rows)


class Results(object):
    """All output files found in a ``SRIM Outputs`` directory.  Missing files
    give ``None`` attributes."""

    _files = (('ioniz', Ioniz), ('phonons', Phonons), ('e2recoil', EnergyToRecoils),
              ('vacancy', Vacancy), ('novac', NoVacancy), ('range', Range),
              ('lateral', Lateral), ('transmit', Transmit), ('backscat', Backscat),
              ('sputter', Sputter), ('range3d', Range3D), ('trimout', TrimOut), ('exyz', EXYZ),
              ('collisions', Collisions))

    def __init__(self, directory):
        self.directory = directory
        for attr, cls in self._files:
            p = os.path.join(directory, cls.name)
            setattr(self, attr, cls(directory) if os.path.exists(p) else None)
