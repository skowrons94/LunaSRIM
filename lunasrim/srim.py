"""pysrim compatible interface to the Python TRIM engine.

The classes mirror those of `pysrim <https://github.com/costrouc/pysrim>`_
(``Ion``, ``Layer``, ``Target``, ``TRIM``) so that existing scripts only
need to change their import::

    from lunasrim.srim import Ion, Layer, Target, TRIM

    ion = Ion('He', energy=2e6)                       # energy in eV
    layer = Layer({'Au': {'stoich': 1.0}}, density=19.3, width=3000.0)
    trim = TRIM(Target([layer]), ion, number_ions=100, calculation=1)
    results = trim.run('output')                     # no SRIM.exe needed

Instead of running ``TRIM.exe`` under wine, :meth:`TRIM.run` executes the
Python re-implementation, which reproduces SRIM-2013 TRIM ion for ion.
"""
import os
import random
import shutil
import warnings

from . import data_files
from .trimin import TrimInput, TrimInputElement, TrimInputLayer
from .trim import TrimEngine
from .output import TrimOutput, TrimDistributions
from .results import Results


def _num(v):
    """number as written into TRIM.IN (integral values without a decimal
    point, as the TRIM user interface writes them)"""
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, int):
        return str(v)
    v = float(v)
    if v.is_integer() and abs(v) < 1e15:
        return str(int(v))
    s = repr(v)
    return s


class Element(object):
    """A chemical element (symbol, atomic number and mass in amu)"""

    def __init__(self, identifier, mass=None, _ion=False):
        if isinstance(identifier, Element):
            z, symbol = identifier.atomic_number, identifier.symbol
            mass = identifier.mass if mass is None else mass
        elif isinstance(identifier, int):
            z = identifier
            symbol = data_files.element(z).symbol
        else:
            z = data_files.symbol_to_z(identifier)
            symbol = data_files.element(z).symbol
        a = data_files.atomdata(z)
        self.atomic_number = z
        self.symbol = symbol
        self.name = a.name
        if mass is None:
            # TRIM's defaults: most abundant isotope for the ion, natural
            # atomic weight for target atoms
            mass = a.mai_weight if _ion else a.natural_weight
        self.mass = mass

    def __repr__(self):
        return '<Element %s Z=%d mass=%s>' % (self.symbol, self.atomic_number, self.mass)

    def __eq__(self, other):
        return isinstance(other, Element) and (self.atomic_number, self.mass) == (other.atomic_number, other.mass)

    def __hash__(self):
        return hash((self.atomic_number, self.mass))


class Ion(Element):
    """The projectile.  ``energy`` is in eV (as in pysrim)."""

    def __init__(self, identifier, energy, mass=None):
        Element.__init__(self, identifier, mass, _ion=True)
        self.energy = energy

    def __repr__(self):
        return '<Ion %s Z=%d mass=%s energy=%s eV>' % (self.symbol, self.atomic_number, self.mass, self.energy)


class Layer(object):
    """A target layer.

    ``elements`` maps an element (symbol, atomic number or :class:`Element`)
    to a dict with ``stoich`` and optionally ``E_d`` (displacement energy),
    ``lattice`` (lattice binding energy) and ``surface`` (surface binding
    energy), all in eV.  ``density`` is in g/cm3 and ``width`` in Angstrom;
    ``phase`` is 0 for a solid and 1 for a gas.
    """

    def __init__(self, elements, density, width, phase=0, name=None, bragg_correction=1.0):
        self.elements = []
        for key, props in elements.items():
            el = Element(key)
            if not isinstance(props, dict):
                props = {'stoich': props}
            p = {'stoich': 1.0, 'E_d': 25.0, 'lattice': 3.0,
                 'surface': data_files.atomdata(el.atomic_number).heat_sublimation}
            p.update(props)
            self.elements.append((el, p))
        self.density = density
        self.width = width
        self.phase = phase
        self.bragg_correction = bragg_correction
        self.name = name or ''.join(el.symbol for el, p in self.elements)

    @classmethod
    def from_formula(cls, formula, density, width, phase=0, name=None):
        """Layer from a chemical formula such as 'SiO2' or 'Ta1N1'"""
        import re
        elements = {}
        for sym, n in re.findall(r'([A-Z][a-z]?)([0-9.]*)', formula):
            elements[sym] = {'stoich': float(n) if n else 1.0}
        return cls(elements, density, width, phase=phase, name=name or formula)


class Target(object):
    """A stack of layers (the first one is hit by the ion)"""

    def __init__(self, layers):
        self.layers = list(layers)

    @property
    def width(self):
        return sum(l.width for l in self.layers)


class TRIM(object):
    """A TRIM calculation (same arguments as pysrim's ``srim.TRIM``).

    ``calculation``: 1 = ion distribution and quick damage (Kinchin-Pease),
    2 = full damage cascades, 3 = monolayer collision steps / sputtering.
    """

    def __init__(self, target, ion, number_ions=1000, calculation=1, angle_ions=0.0,
                 random_seed=None, reminders=0, autosave=10000, plot_mode=5,
                 plot_xmin=0.0, plot_xmax=0.0, ranges=False, backscattered=False,
                 transmit=False, sputtered=False, collisions=False, exyz=0,
                 bragg_correction=1.0, description='lunasrim calculation', **kwargs):
        if 'collissions' in kwargs:          # spelling used in some pysrim scripts
            collisions = kwargs.pop('collissions')
        if kwargs:
            raise TypeError('unknown arguments: %s' % ', '.join(sorted(kwargs)))
        self.target = target
        self.ion = ion
        self.settings = dict(
            number_ions=number_ions, calculation=calculation, angle_ions=angle_ions,
            random_seed=random.randint(0, 100000) if random_seed is None else random_seed,
            reminders=reminders, autosave=autosave, plot_mode=plot_mode,
            plot_xmin=plot_xmin, plot_xmax=plot_xmax, ranges=ranges,
            backscattered=backscattered, transmit=transmit, sputtered=sputtered,
            collisions=collisions, exyz=exyz, bragg_correction=bragg_correction,
            description=description)

    # ------------------------------------------------------------- TRIM.IN
    def trim_input_text(self):
        """the TRIM.IN file of this calculation"""
        s = self.settings
        ion, layers = self.ion, self.target.layers
        elements = []
        for l in layers:
            for el, p in l.elements:
                if el not in [e for e, q in elements]:
                    elements.append((el, p))
        ti = TrimInput()
        ti.title = s['description']
        ti.ion_z = ion.atomic_number
        ti.ion_mass = ion.mass
        ti.energy = ion.energy / 1000.0
        ti.angle = s['angle_ions']
        ti.number_ions = s['number_ions']
        ti.bragg_ion = s['bragg_correction']
        ti.autosave = s['autosave']
        ti.cascade = s['calculation']
        ti.seed1 = s['random_seed']
        ti.seed2 = s['reminders']
        ti.file_ranges = int(bool(s['ranges']))
        ti.file_backscatt = int(bool(s['backscattered']))
        ti.file_transmit = int(bool(s['transmit']))
        ti.file_sputtered = int(bool(s['sputtered']))
        ti.file_collisions = int(s['collisions'])
        ti.file_exyz = s['exyz']
        ti.target_name = ''.join(l.name for l in layers)[:30] or 'target'
        ti.plot_type = s['plot_mode']
        ti.xmin = s['plot_xmin']
        ti.xmax = s['plot_xmax']
        for el, p in elements:
            ti.elements.append(TrimInputElement(el.symbol, el.atomic_number, el.mass))
            ti.ed.append(p['E_d'])
            ti.lattice.append(p['lattice'])
            ti.surface.append(p['surface'])
        for l in layers:
            stoich = []
            for el, p in elements:
                st = 0
                for e2, p2 in l.elements:
                    if e2 == el:
                        st = p2['stoich']
                stoich.append(st)
            tl = TrimInputLayer(l.name, l.width, l.density, stoich)
            tl.phase = l.phase
            tl.bragg = l.bragg_correction
            ti.layers.append(tl)
        return _trim_in_text(ti)

    def trim_input(self):
        """:class:`TrimInput` as TRIM would read the generated TRIM.IN"""
        return TrimInput.parse(self.trim_input_text())

    # ----------------------------------------------------------------- run
    def run(self, directory='.', progress=None):
        """Run the calculation with the Python engine.

        TRIM.IN is written into ``directory`` and the output files into
        ``directory/SRIM Outputs`` (the layout of a SRIM installation).
        Returns a :class:`lunasrim.results.Results`.
        """
        if self.settings['collisions']:
            warnings.warn('COLLISON.txt is not produced by lunasrim; the setting is ignored')
        if self.settings['exyz']:
            warnings.warn('EXYZ.txt is not produced by lunasrim; the setting is ignored')
        os.makedirs(directory, exist_ok=True)
        text = self.trim_input_text()
        with open(os.path.join(directory, 'TRIM.IN'), 'wb') as fh:
            fh.write(text.encode('latin-1'))
        outdir = os.path.join(directory, 'SRIM Outputs')
        run_trim(TrimInput.parse(text), outdir, progress=progress)
        return Results(outdir)

    @staticmethod
    def copy_output_files(src_directory, dest_directory, check_srim_output=True):
        """Copy TRIM.IN and the output files of a run (pysrim compatible)"""
        known = ('TRIM.IN', 'PHONON.txt', 'E2RECOIL.txt', 'IONIZ.txt', 'LATERAL.txt',
                 'NOVAC.txt', 'RANGE.txt', 'VACANCY.txt', 'COLLISON.txt', 'BACKSCAT.txt',
                 'SPUTTER.txt', 'RANGE_3D.txt', 'TRANSMIT.txt', 'TRIMOUT.txt', 'TDATA.txt')
        os.makedirs(dest_directory, exist_ok=True)
        found = False
        for d in (src_directory, os.path.join(src_directory, 'SRIM Outputs')):
            if not os.path.isdir(d):
                continue
            for name in os.listdir(d):
                if name in known:
                    shutil.copy(os.path.join(d, name), os.path.join(dest_directory, name))
                    found = found or name.endswith('.txt')
        if check_srim_output and not found:
            raise ValueError('no SRIM output files found in %s' % src_directory)


def _trim_in_text(ti):
    """TRIM.IN text with TRIM's number formatting"""
    nl = '\r\n'
    o = []
    o.append(ti.title)
    o.append('Ion: Z, Mass [amu], Energy [keV], Angle [degrees], Number Ions, Bragg Corr, AutoSave Number')
    o.append(' '.join(_num(v) for v in (ti.ion_z, ti.ion_mass, ti.energy, ti.angle,
                                        ti.number_ions, ti.bragg_ion, ti.autosave)))
    o.append('Cascades(1=Kitchn-Peese, 2=Full-Cascade, 3=Sputtering, 4-5=Ions;6-7=Neutrons), '
             'Random Number Seed, Reminders')
    o.append('%d %d %d' % (ti.cascade, ti.seed1, ti.seed2))
    o.append('Diskfiles (0=no,1=yes): RANGES.txt, BACKSCATT.txt, TRANSMIT.txt, Sputtered, '
             'COLLISIONS.txt(0=no, 1=Ion, 2=Ion+Recoils), Special EXYZ.txt file')
    o.append(' '.join(_num(v) for v in (ti.file_ranges, ti.file_backscatt, ti.file_transmit,
                                        ti.file_sputtered, ti.file_collisions, ti.file_exyz)))
    o.append('Target material : Number of Elements, Number of Layers')
    o.append('"%s" %d %d' % (ti.target_name, len(ti.elements), len(ti.layers)))
    o.append('PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]')
    o.append('%d %s %s' % (ti.plot_type, _num(ti.xmin), _num(ti.xmax)))
    o.append('Target Elements:    Z   Mass [amu]')
    for i, e in enumerate(ti.elements):
        o.append('Atom %d = %s =     %d %s' % (i + 1, e.symbol, e.z, _num(e.mass)))
    o.append('Layer    Layer Name   Width Density  ' + ' '.join('%s(%d)' % (e.symbol, e.z) for e in ti.elements))
    o.append('Number   Description  (Ang) (g/cm^3)' + ' '.join('Stoich' for e in ti.elements))
    for j, l in enumerate(ti.layers):
        o.append('%d "%s" %s %s %s' % (j + 1, l.name, _num(l.width), _num(l.density),
                                       ' '.join(_num(s) for s in l.stoich)))
    o.append('0  Target layer phases (0=Solid, 1=Gas)')
    o.append(' '.join('%d' % l.phase for l in ti.layers))
    o.append('Target Compound Corrections (Bragg)')
    o.append(' ' + ' '.join(_num(l.bragg) for l in ti.layers))
    o.append('Individual target atom displacement energies (eV)')
    o.append(' ' + ' '.join(_num(x) for x in ti.ed))
    o.append('Individual target atom lattice binding energies (eV)')
    o.append(' ' + ' '.join(_num(x) for x in ti.lattice))
    o.append('Individual target atom surface binding energies (eV)')
    o.append(' ' + ' '.join(_num(x) for x in ti.surface))
    o.append('Stopping Power Version (1=2011, 0=2011)')
    o.append('%d' % ti.stopping_version)
    return nl.join(o) + nl


def run_trim(ti, directory, number_ions=None, progress=None, date=None):
    """Run a TRIM calculation described by a :class:`TrimInput` and write
    all output files into ``directory``.  Returns the :class:`TrimEngine`.

    ``progress`` may be a callable ``progress(ion_number, engine)`` called
    after every ion.
    """
    eng = TrimEngine(ti)
    eng.out = TrimOutput(eng, directory, date=date)
    eng.out.start()
    n = int(ti.number_ions if number_ions is None else number_ions)
    try:
        for i in range(n):
            eng.run_ion()
            if progress is not None:
                progress(i + 1, eng)
    finally:
        eng.out.close()
    TrimDistributions(eng, directory).write_all()
    return eng
