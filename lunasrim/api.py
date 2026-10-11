"""Simple Python API: ions and targets by name, energies in keV.

The functions of this module are re-exported by the package::

    import lunasrim

    rows = lunasrim.stopping_table('H', 'Si', emin=10, emax=10000)
    rows = lunasrim.stopping_table('He', 'SiO2', density=2.32, emin=100, emax=5000, units='eV/A')
    r = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um')
    r = lunasrim.energy_loss('He', 'N2', energy=300, thickness=2e19, unit='atoms/cm2', gas=True)
    res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=1000)
    m = lunasrim.material('Kapton')

They only resolve names and defaults; the calculations are done by the
low-level functions (:func:`lunasrim.sr.stopping_table`,
:func:`lunasrim.straggling.energy_loss`, :class:`lunasrim.srim.TRIM`), so
the results are identical to theirs.

**Ions** are given by element symbol (``'He'``), element name
(``'Helium'``) or atomic number (``2``).  The mass defaults to the most
abundant isotope (SRIM's default, e.g. 4.003 amu for He); ``mass=`` sets
another one.

**Targets** (see :func:`material`) are an element (``'Si'``), a chemical
formula (``'SiO2'``, ``'Al2O3'``, ``'C8H8'``), a name from SRIM's compound
dictionary (``'Kapton'``, ``'Mylar'``, ``'Water_Liquid'``), a dict
``{'Si': 1, 'O': 2}`` or a :class:`Material`.
"""
import math
import os
import re
import tempfile

from . import compounds, data_files, sr, straggling
from .stopping import GAS_ELEMENTS

__all__ = ['Material', 'material', 'resolve_ion', 'stopping_table', 'energy_loss', 'trim',
           'stopping_units', 'parse_energy', 'parse_thickness', 'STOPPING_UNITS']

# readable names of SRIM's stopping units (sr.UNITS)
STOPPING_UNITS = {1: 'eV/A', 2: 'keV/um', 3: 'MeV/mm', 4: 'keV/(ug/cm2)', 5: 'MeV/(mg/cm2)',
                  6: 'keV/(mg/cm2)', 7: 'eV/(1E15 atoms/cm2)', 8: 'LSS reduced units'}

_UNIT_ALIASES = {
    'ev/a': 1, 'ev/angstrom': 1,
    'kev/um': 2, 'kev/micron': 2,
    'mev/mm': 3,
    'kev/ug/cm2': 4,
    'mev/mg/cm2': 5,
    'kev/mg/cm2': 6,
    'ev/1e15atoms/cm2': 7, 'ev/1e15/cm2': 7, 'ev/1e15at/cm2': 7,
    'lss': 8, 'l.s.s.': 8, 'lssreducedunits': 8, 'l.s.s.reducedunits': 8, 'reduced': 8,
}

# default displacement energies (eV) of the TRIM user interface (as lunasrim.web.core)
_DISPLACEMENT = {1: 10, 6: 28, 7: 28, 8: 28, 14: 15}
_LATTICE = 3.0


# ------------------------------------------------------------------ helpers
def _elements_by_symbol():
    return {e.symbol: z for z, e in data_files.load_atomdata().items()}


def _z_of(name):
    """atomic number of a symbol, name or number (int or digit string);
    None if it is not an element"""
    if isinstance(name, bool):
        return None
    if isinstance(name, int):
        return name if 1 <= name <= 92 else None
    s = str(name).strip()
    if s.isdigit():
        z = int(s)
        return z if 1 <= z <= 92 else None
    atoms = data_files.load_atomdata()
    for z, a in atoms.items():
        if a.symbol == s:
            return z
    low = s.lower()
    for z, a in atoms.items():
        if a.symbol.lower() == low or a.name.lower() == low:
            return z
    return None


def resolve_ion(ion, mass=None):
    """``(Z, mass)`` of an ion given by symbol, name or atomic number.

    The mass defaults to the most abundant isotope (amu), as in SRIM and in
    the web application.  Raises ``ValueError`` for an unknown element.
    """
    z = _z_of(ion)
    if z is None:
        raise ValueError('unknown ion %r (give an element symbol, name or Z = 1-92)' % (ion,))
    if mass is None:
        mass = data_files.atomdata(z).mai_weight
    elif not mass > 0:
        raise ValueError('the ion mass must be positive')
    return z, float(mass)


class Material(object):
    """A target material resolved to the numbers SRIM needs.

    Attributes:

    * ``name``: the name given (or of the compound found);
    * ``elements``: list of ``(Z, stoichiometry, mass)``, mass in amu
      (natural atomic weight by default);
    * ``density``: g/cm3;
    * ``gas``: True for a gas target (SRIM's gas phase stopping);
    * ``correction``: compound (Bragg) correction factor, 1 by default;
    * ``source``: ``'element'``, ``'formula'``, ``'compound'`` or ``'custom'``.
    """

    def __init__(self, name, elements, density, gas=False, correction=1.0, source='custom'):
        self.name = name
        # integral stoichiometries stay integers: as in SR.IN / TRIM.IN files,
        # SRIM's Visual Basic arithmetic then works in Single precision
        self.elements = [(int(z), _stoich(st), float(m)) for (z, st, m) in elements]
        self.density = density
        self.gas = bool(gas)
        self.correction = float(correction)
        self.source = source

    @property
    def phase(self):
        """0 for a solid, 1 for a gas (the SRIM input value)"""
        return 1 if self.gas else 0

    @property
    def symbols(self):
        return [data_files.atomdata(z).symbol for z, st, m in self.elements]

    @property
    def mean_mass(self):
        """stoichiometric mean atomic mass (amu)"""
        tot = sum(st for z, st, m in self.elements)
        return sum(st / tot * m for z, st, m in self.elements)

    @property
    def atom_density(self):
        """atoms/cm3"""
        return 6.022e23 * self.density / self.mean_mass

    def formula(self):
        """composition as a formula, e.g. 'Si1 O2'"""
        return ' '.join('%s%s' % (s, ('%.6g' % st)) for s, (z, st, m) in zip(self.symbols, self.elements))

    def as_dict(self):
        return {'name': self.name, 'density': self.density, 'gas': self.gas,
                'correction': self.correction, 'source': self.source,
                'elements': [{'z': z, 'symbol': s, 'stoich': st, 'mass': m}
                             for s, (z, st, m) in zip(self.symbols, self.elements)]}

    def __repr__(self):
        return '<Material %r: %s, %.6g g/cm3, %s>' % (self.name, self.formula(), self.density,
                                                    'gas' if self.gas else 'solid')


def _stoich(v):
    if isinstance(v, bool):
        raise ValueError('not a valid stoichiometry: %r' % (v,))
    if isinstance(v, int):
        return v
    v = float(v)
    if not (v >= 0 and math.isfinite(v)):
        raise ValueError('not a valid stoichiometry: %r' % (v,))
    return v


_FORMULA = re.compile(r'([A-Z][a-z]?)(\d*\.?\d*)')


def _parse_formula(text):
    """[(Z, n)] of a formula such as 'SiO2' or 'C8H8', None if it is not one"""
    s = text.replace(' ', '')
    if not s or not s[0].isupper():
        return None
    pos, out = 0, []
    for m in _FORMULA.finditer(s):
        if m.start() != pos:
            return None
        pos = m.end()
        sym, n = m.group(1), m.group(2)
        z = _elements_by_symbol().get(sym)
        if z is None:
            return None
        try:
            n = (float(n) if '.' in n else int(n)) if n else 1
        except ValueError:
            return None
        if n <= 0:
            return None
        out.append((z, n))
    if pos != len(s) or not out:
        return None
    merged = []
    for z, n in out:                 # 'CH3COOH' -> C2 H4 O2
        for i, (z2, n2) in enumerate(merged):
            if z2 == z:
                merged[i] = (z, n2 + n)
                break
        else:
            merged.append((z, n))
    return merged


def _norm(s):
    return re.sub(r'[\s_\-,]+', ' ', s.lower()).strip()


def _strip_ref(s):
    """'Water_Liquid (ICRU-276)' -> 'Water_Liquid'"""
    return re.sub(r'\s*\([^)]*\)\s*$', '', s).strip()


def _same(c1, c2):
    return c1.density == c2.density and c1.elements == c2.elements and c1.gas == c2.gas


def _find_compound(name, fuzzy):
    """compound of SRIM's dictionary called ``name`` (case-insensitive).
    Exact match of the name or of the name without its '(ICRU-...)'
    reference; with ``fuzzy`` also names or titles that start with it."""
    q = _norm(name)
    cs = compounds.load()
    for key in (lambda c: _norm(c.name), lambda c: _norm(_strip_ref(c.name))):
        for c in cs:
            if key(c) == q:
                return c
    if not fuzzy:
        return None
    found = []
    for words in (False, True):
        # names or descriptions starting with the text, then containing it as words
        for c in cs:
            for text in (c.name, c.title):
                t = _norm(text)
                hit = (' %s ' % q) in (' %s ' % re.sub(r'[()]', ' ', t)) if words else \
                    (t == q or t.startswith(q + ' ') or t.startswith(q + '('))
                if hit:
                    if not any(_same(c, f) or c.name == f.name for f in found):
                        found.append(c)
                    break
        if found:
            break
    if len(found) == 1:
        return found[0]
    if len(found) > 1:
        raise ValueError('%r matches several compounds of the SRIM dictionary: %s; give the full name'
                         % (name, '; '.join('%r (%g g/cm3)' % (c.name, c.density) for c in found[:12])))
    return None


def _compounds_like(els):
    """compounds with the composition of the formula ``els`` [(Z, n)]"""
    tot = sum(n for z, n in els)
    want = sorted((z, n / tot) for z, n in els)
    out = []
    for c in compounds.load():
        have = sorted(c.elements)
        if len(have) == len(want) and all(z1 == z2 and abs(f1 - f2) < 2e-3
                                          for (z1, f1), (z2, f2) in zip(have, want)):
            if c.name not in [o.name for o in out]:
                out.append(c)
    return out


def _element_material(z, name, stoich, density, gas):
    a = data_files.atomdata(z)
    sc = data_files.element(z)
    natural_gas = z in GAS_ELEMENTS
    gas = natural_gas if gas is None else bool(gas)
    if density is None:
        if gas:
            if not sc.gas_density > 0:
                raise ValueError('SRIM has no gas density for %s; give density= (g/cm3)' % a.symbol)
            density = sc.gas_density
        else:
            density = a.density
    return Material(name, [(z, stoich, a.natural_weight)], density, gas, source='element')


def material(target, density=None, gas=None, correction=None):
    """Resolve a target to a :class:`Material`.

    ``target`` is one of (tried in this order):

    1. an element symbol (``'Si'``, case-sensitive first, so ``'CO'`` is
       carbon monoxide and not cobalt), element name (``'silicon'``) or
       atomic number;
    2. the name of a compound of SRIM's dictionary (``Compound.dat``,
       356 entries), case-insensitive, with or without its reference in
       parentheses: ``'Mylar'``, ``'Water_Liquid'``, ``'Kapton Polyimide
       Film (ICRU-179)'``;
    3. a chemical formula: ``'SiO2'``, ``'Al2O3'``, ``'C8H8'``, ``'H2O'``
       (a formula of one element, such as ``'N2'``, is that element);
    4. the beginning of a compound name or description, or else words of
       it, if they identify a single entry: ``'Kapton'``, ``'Teflon'``,
       ``'quartz'`` (an ambiguous text such as ``'Water'`` gives a
       ``ValueError`` listing the candidates);
    5. a dict ``{'Si': 1, 'O': 2}`` (symbol or Z to stoichiometry) or a
       :class:`Material` (returned with the overrides applied).

    Defaults:

    * ``density`` (g/cm3): for an element the solid density of SRIM's
      tables, or SRIM's gas density for the gaseous elements (H, He, N, O,
      F, Ne, Cl, Ar, Br, Kr, I, Xe, At, Rn), as the web application does;
      for a compound its dictionary density.  It is required for a formula
      of several elements and for a dict (``ValueError`` otherwise).
    * ``gas``: True for gaseous elements and gas compounds, False
      otherwise.  It selects SRIM's gas phase stopping.
    * ``correction``: compound (Bragg) correction, 1.

    Atomic masses are the natural atomic weights.
    """
    if isinstance(target, Material):
        m = Material(target.name, target.elements, target.density, target.gas, target.correction,
                     target.source)
        if density is not None:
            m.density = density
        if gas is not None:
            m.gas = bool(gas)
    elif isinstance(target, dict):
        els = []
        for k, st in target.items():
            z = _z_of(k)
            if z is None:
                raise ValueError('unknown element %r' % (k,))
            els.append((z, _stoich(st), data_files.atomdata(z).natural_weight))
        if not els:
            raise ValueError('the target has no elements')
        if len(els) == 1 and density is None:
            m = _element_material(els[0][0], ''.join(str(k) for k in target), els[0][1], None, gas)
        else:
            if density is None:
                raise ValueError('give the density (g/cm3) of the target %r' % (target,))
            m = Material(''.join('%s%g' % (data_files.atomdata(z).symbol, st) for z, st, ms in els),
                         els, density, bool(gas), source='custom')
    else:
        text = str(target).strip()
        if not text:
            raise ValueError('no target given')
        m = None
        z = None
        if isinstance(target, int) or text in _elements_by_symbol() or text.isdigit() or \
                text.lower() in [a.name.lower() for a in data_files.load_atomdata().values()]:
            z = _z_of(target)
        if z is not None:
            m = _element_material(z, data_files.atomdata(z).symbol, 1, density, gas)
        if m is None:
            c = _find_compound(text, fuzzy=False)
            if c is not None:
                m = _compound_material(c, density, gas)
        if m is None:
            els = _parse_formula(text)
            if els is not None:
                if len(els) == 1:
                    m = _element_material(els[0][0], text, els[0][1], density, gas)
                else:
                    if density is None:
                        like = _compounds_like(els)
                        hint = ''
                        if like:
                            hint = ' (compounds with this composition: %s)' % '; '.join(
                                '%r %g g/cm3' % (c.name, c.density) for c in like[:6])
                        raise ValueError('give the density (g/cm3) of %s%s' % (text, hint))
                    m = Material(text, [(z, n, data_files.atomdata(z).natural_weight) for z, n in els],
                                 density, bool(gas), source='formula')
        if m is None:
            z = _z_of(text)
            if z is not None:
                m = _element_material(z, data_files.atomdata(z).symbol, 1, density, gas)
        if m is None:
            c = _find_compound(text, fuzzy=True)
            if c is not None:
                m = _compound_material(c, density, gas)
        if m is None:
            raise ValueError('unknown target %r: not an element, a chemical formula (e.g. SiO2) '
                             'or a compound of the SRIM dictionary' % (target,))
    if correction is not None:
        m.correction = float(correction)
    if not (m.density > 0):
        raise ValueError('the density must be positive')
    return m


def _compound_material(c, density, gas):
    atoms = data_files.load_atomdata()
    els = [(z, f, atoms[z].natural_weight) for z, f in c.elements]
    return Material(c.name, els, c.density if density is None else density,
                    c.gas if gas is None else bool(gas), source='compound')


# -------------------------------------------------------------------- units
def stopping_units(units):
    """SRIM stopping unit number (1-8) from a number or a name.

    Names (case and spaces do not matter): ``'eV/A'`` (1), ``'keV/um'``
    (2), ``'MeV/mm'`` (3), ``'keV/(ug/cm2)'`` (4), ``'MeV/(mg/cm2)'`` (5),
    ``'keV/(mg/cm2)'`` (6), ``'eV/(1E15 atoms/cm2)'`` (7), ``'LSS'``
    (8, LSS reduced units).  'Angstrom', 'micron' and 'µ' are accepted.
    """
    if isinstance(units, bool):
        raise ValueError('unknown stopping units %r' % (units,))
    if isinstance(units, int) or (isinstance(units, float) and units == int(units)):
        u = int(units)
    else:
        s = str(units).strip().lower().replace('µ', 'u').replace('μ', 'u')
        s = s.replace('Å', 'a').replace('å', 'a').replace('angstroms', 'a').replace('angstrom', 'a')
        s = s.replace('microns', 'um').replace('micron', 'um')
        s = re.sub(r'[\s()]', '', s)
        if s.isdigit():
            u = int(s)
        elif s in _UNIT_ALIASES:
            u = _UNIT_ALIASES[s]
        else:
            raise ValueError('unknown stopping units %r; use 1-8 or one of: %s'
                             % (units, ', '.join(STOPPING_UNITS.values())))
    if u not in STOPPING_UNITS:
        raise ValueError('stopping units must be 1-8')
    return u


_NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
_ENERGY_UNITS = {'ev': 1e-3, 'kev': 1.0, 'mev': 1e3, 'gev': 1e6}


def parse_energy(text):
    """Energy in keV from a number (keV) or a string such as ``'400'``,
    ``'400keV'``, ``'2 MeV'``, ``'50eV'``, ``'1GeV'``."""
    if isinstance(text, (int, float)) and not isinstance(text, bool):
        return float(text)
    m = re.match(r'^\s*(%s)\s*([a-zA-Z]*)\s*$' % _NUMBER, str(text))
    if not m:
        raise ValueError('not a valid energy: %r (examples: 400, 400keV, 2MeV)' % (text,))
    unit = m.group(2).lower() or 'kev'
    if unit not in _ENERGY_UNITS:
        raise ValueError('unknown energy unit %r (eV, keV, MeV, GeV)' % m.group(2))
    return float(m.group(1)) * _ENERGY_UNITS[unit]


# thickness unit aliases -> (unit of straggling.THICKNESS_UNITS, factor)
_THICKNESS_ALIASES = {
    'a': ('A', 1), 'ang': ('A', 1), 'angstrom': ('A', 1), 'angstroms': ('A', 1),
    'nm': ('nm', 1),
    'um': ('um', 1), 'micron': ('um', 1), 'microns': ('um', 1),
    'mm': ('mm', 1), 'cm': ('mm', 10),
    'ug/cm2': ('ug/cm2', 1), 'mg/cm2': ('mg/cm2', 1),
    '1e15/cm2': ('1e15/cm2', 1), '1e15atoms/cm2': ('1e15/cm2', 1), '1e15at/cm2': ('1e15/cm2', 1),
    'atoms/cm2': ('atoms/cm2', 1), 'at/cm2': ('atoms/cm2', 1), '/cm2': ('atoms/cm2', 1),
}


def thickness_unit(unit):
    """``(unit, factor)``: the :data:`lunasrim.straggling.THICKNESS_UNITS`
    name of a thickness unit and the factor to apply to the value (only 'cm'
    has one: 1 cm = 10 mm)"""
    s = str(unit).strip().replace('µ', 'u').replace('μ', 'u')
    s = s.replace('Å', 'A').replace('å', 'A')
    s = re.sub(r'\s+', '', s)
    low = s.lower().replace('^2', '2')
    if low in _THICKNESS_ALIASES:
        return _THICKNESS_ALIASES[low]
    raise ValueError('unknown thickness unit %r; use one of: A, nm, um, mm, cm, ug/cm2, mg/cm2, '
                     '1e15/cm2, atoms/cm2' % (unit,))


def parse_thickness(text, unit=None):
    """``(value, unit)`` of a thickness for :func:`energy_loss` from a
    string with a unit, e.g. ``'1um'``, ``'200 nm'``, ``'50A'``,
    ``'100ug/cm2'``, ``'0.5mg/cm2'``, ``'2e19atoms/cm2'``, ``'20e15/cm2'``,
    ``'22cm'``.  A plain number takes ``unit`` (required then).  The unit
    returned is one of :data:`lunasrim.straggling.THICKNESS_UNITS`."""
    if isinstance(text, (int, float)) and not isinstance(text, bool):
        value, u = float(text), unit
    else:
        m = re.match(r'^\s*(%s)\s*(.*?)\s*$' % _NUMBER, str(text))
        if not m:
            raise ValueError('not a valid thickness: %r (examples: 1um, 200nm, 100ug/cm2, 2e19atoms/cm2)'
                             % (text,))
        value, u = float(m.group(1)), m.group(2) or unit
    if not u:
        raise ValueError('the thickness %r needs a unit (A, nm, um, mm, cm, ug/cm2, mg/cm2, 1e15/cm2, '
                         'atoms/cm2)' % (text,))
    name, factor = thickness_unit(u)
    return (value * factor if factor != 1 else value), name


# -------------------------------------------------------------- calculations
def _sr_input(ion, target, mass, density, gas, correction, units, emin, emax):
    z1, m1 = resolve_ion(ion, mass)
    mat = material(target, density, gas, correction)
    si = sr.SrInput()
    si.ion_z, si.ion_mass = z1, m1
    si.phase, si.density, si.correction = mat.phase, mat.density, mat.correction
    si.elements = [sr.SrElement(z, data_files.atomdata(z).name, st, m) for (z, st, m) in mat.elements]
    si.units = stopping_units(units)
    si.emin, si.emax = parse_energy(emin), parse_energy(emax)
    si.output_name = ('%s in %s' % (data_files.atomdata(z1).name, mat.name))[:60]
    return si, mat


def _check_energies(emin, emax):
    if not (0 < emin <= emax):
        raise ValueError('need 0 < emin <= emax (keV)')


def stopping_table(ion, target, emin=10.0, emax=10000.0, units='MeV/(mg/cm2)', density=None,
                   gas=None, mass=None, correction=None):
    """Stopping and range table (SRIM's "Stopping and Range Tables", SRModule.exe).

    ``ion``: symbol, name or Z (``mass`` in amu, default most abundant
    isotope); ``target``: element, formula or compound (see
    :func:`material`, which also gives the defaults of ``density`` g/cm3,
    ``gas`` and ``correction``).  ``emin`` and ``emax`` are in keV (a
    string such as ``'2MeV'`` also works); the table uses SRIM's energy
    grid between them.  ``units``: stopping units, 1-8 or a name such as
    ``'eV/A'``, ``'keV/um'``, ``'MeV/(mg/cm2)'`` (default, SRIM's default),
    ``'eV/(1E15 atoms/cm2)'``, see :func:`stopping_units`.

    Returns a list of tuples, one per energy::

        (energy keV, S_electronic, S_nuclear, projected range A,
         longitudinal straggling A, lateral straggling A)

    with the stopping powers in ``units`` and the ranges in Angstrom.  It
    is the same list as :func:`lunasrim.sr.stopping_table` gives (which is
    identical to SRModule.exe's table).  Raises ``ValueError`` for bad
    input and :class:`lunasrim.sr.SrError` where SRModule.exe fails.
    """
    si, mat = _sr_input(ion, target, mass, density, gas, correction, units, emin, emax)
    _check_energies(si.emin, si.emax)
    return sr.stopping_table(si.ion_z, si.ion_mass, mat.elements, mat.density, si.emin, si.emax,
                             phase=mat.phase, correction=mat.correction, units=si.units)


def stopping_report(ion, target, emin=10.0, emax=10000.0, units='MeV/(mg/cm2)', density=None,
                    gas=None, mass=None, correction=None):
    """The calculated :class:`lunasrim.sr.StoppingTable` (same arguments as
    :func:`stopping_table`); ``table.lines()`` is the SRModule.exe output
    file and ``table.si`` the SR.IN input."""
    si, mat = _sr_input(ion, target, mass, density, gas, correction, units, emin, emax)
    _check_energies(si.emin, si.emax)
    t = sr.StoppingTable(si, strict=False)
    t.calculate()
    return t


def energy_loss(ion, target, energy, thickness, unit='A', angle=0.0, density=None, gas=None,
                mass=None, correction=None, steps=500):
    """Energy loss and energy straggling of a beam crossing one layer.

    ``energy`` in keV (or a string such as ``'2MeV'``); ``thickness`` in
    ``unit``: 'A', 'nm', 'um', 'mm', 'cm', 'ug/cm2', 'mg/cm2', '1e15/cm2'
    or 'atoms/cm2' (a string with its unit such as ``'1um'`` also works);
    ``angle``: tilt of the beam from the layer normal in degrees (the path
    length is thickness / cos(angle)).  ``ion``, ``target``, ``density``,
    ``gas``, ``mass`` and ``correction`` as in :func:`stopping_table`.

    Returns the dict of :func:`lunasrim.straggling.energy_loss` (same
    numbers), with in particular:

    * ``energy_out``: mean exit energy (keV); ``energy_loss``,
      ``electronic_loss``, ``nuclear_loss`` (keV); ``mean_energy`` in the
      layer;
    * ``sigma``, ``fwhm``: energy straggling (keV), TRIM.exe's model;
      ``sigma_bohr``, ``fwhm_bohr``: Bohr's formula;
    * ``thickness`` (path length, A), ``areal_density`` (atoms/cm2),
      ``mass_thickness`` (ug/cm2), ``atom_density`` (atoms/cm3);
    * ``stopped``: True if the ions stop in the layer;
    * ``profile``: ``depth`` (A), ``energy``, ``sigma``, ``sigma_bohr``
      (keV), ``se``, ``sn`` (eV/A) along the path;

    plus ``angle`` and ``path_factor`` (1/cos(angle)) and ``material``
    (the :class:`Material` used).
    """
    z1, m1 = resolve_ion(ion, mass)
    mat = material(target, density, gas, correction)
    if isinstance(thickness, str):
        value, u = parse_thickness(thickness, unit)
    else:
        value, u = parse_thickness(float(thickness), unit)
    e = parse_energy(energy)
    if not (-89.99 <= angle <= 89.99):
        raise ValueError('the angle must be within +-89.99 degrees')
    tilt = 1 / math.cos(math.radians(angle))
    r = straggling.energy_loss(z1, m1, e, mat.elements, mat.density, value * tilt, u,
                               phase=mat.phase, correction=mat.correction, steps=steps)
    r['angle'], r['path_factor'], r['material'] = float(angle), tilt, mat
    return r


def trim(ion, target, energy, thickness, unit='A', ions=1000, calculation=1, angle=0.0, seed=None,
         directory=None, density=None, gas=None, mass=None, correction=None, progress=None):
    """Run a TRIM Monte Carlo calculation of an ion beam in one layer.

    ``energy`` in keV; ``thickness`` in ``unit`` (as in
    :func:`energy_loss`, converted to Angstrom with the layer density);
    ``ions``: number of ions; ``calculation``: 1 = ion distribution and
    quick Kinchin-Pease damage, 2 = full damage cascades, 3 = monolayer
    collision steps / surface sputtering; ``angle``: incidence angle
    (degrees); ``seed``: random number seed (random if None; the same seed
    gives the same ions).  ``ion``, ``target``, ``density``, ``gas``,
    ``mass`` and ``correction`` as in :func:`stopping_table`.

    The target atoms get the TRIM user interface defaults, as in the web
    application: displacement energy 10 eV (H), 28 eV (C, N, O), 15 eV
    (Si), 25 eV otherwise; lattice binding energy 3 eV; surface binding
    energy = heat of sublimation (SRIM's ATOMDATA).  RANGE_3D, BACKSCAT,
    TRANSMIT and SPUTTER files are switched on; COLLISON and EXYZ are not
    (use :class:`lunasrim.srim.TRIM` for full control).

    TRIM.IN is written into ``directory`` and the output files into
    ``directory/SRIM Outputs``.  With ``directory=None`` a new temporary
    directory is made (and kept).  ``progress(ion_number, engine)`` is
    called after every ion if given.

    Returns :class:`lunasrim.results.Results` (``.directory`` is the
    ``SRIM Outputs`` path), e.g. ``res.range.average_range`` (A),
    ``res.vacancy.vacancies_per_ion``, ``len(res.transmit)``.
    """
    from .srim import Element, Ion, Layer, Target, TRIM
    z1, m1 = resolve_ion(ion, mass)
    mat = material(target, density, gas, correction)
    value, u = parse_thickness(thickness, unit)
    width = straggling.to_angstrom(value, u, mat.density, mat.mean_mass)
    if not width > 0:
        raise ValueError('the thickness must be positive')
    e = parse_energy(energy)
    if not e > 0:
        raise ValueError('the ion energy must be positive')
    if int(ions) < 1:
        raise ValueError('the number of ions must be at least 1')
    if calculation not in (1, 2, 3):
        raise ValueError('calculation must be 1, 2 or 3')
    if not (-89.99 <= angle <= 89.99):
        raise ValueError('the angle must be within +-89.99 degrees')
    els = {}
    for z, st, m in mat.elements:
        els[Element(z, m)] = {'stoich': st, 'E_d': _DISPLACEMENT.get(z, 25), 'lattice': _LATTICE,
                              'surface': data_files.atomdata(z).heat_sublimation}
    layer = Layer(els, density=mat.density, width=width, phase=mat.phase,
                  name=mat.name[:20].replace('"', ''), bragg_correction=mat.correction)
    calc = TRIM(Target([layer]), Ion(z1, energy=e * 1000.0, mass=m1), number_ions=int(ions),
                calculation=calculation, angle_ions=angle, random_seed=seed, ranges=True,
                backscattered=True, transmit=True, sputtered=True,
                description=('%s (%g keV) into %s' % (data_files.atomdata(z1).symbol, e, mat.name))[:60])
    if directory is None:
        directory = tempfile.mkdtemp(prefix='lunasrim_')
    return calc.run(os.path.abspath(directory), progress=progress)
