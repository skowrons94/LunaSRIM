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


class Results(object):
    """All output files found in a ``SRIM Outputs`` directory.  Missing files
    give ``None`` attributes."""

    _files = (('ioniz', Ioniz), ('phonons', Phonons), ('e2recoil', EnergyToRecoils),
              ('vacancy', Vacancy), ('novac', NoVacancy), ('range', Range),
              ('lateral', Lateral), ('transmit', Transmit), ('backscat', Backscat),
              ('sputter', Sputter), ('range3d', Range3D), ('trimout', TrimOut))

    def __init__(self, directory):
        self.directory = directory
        for attr, cls in self._files:
            p = os.path.join(directory, cls.name)
            setattr(self, attr, cls(directory) if os.path.exists(p) else None)
