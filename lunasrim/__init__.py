"""lunasrim: a Python 3 re-implementation of SRIM-2013 (TRIM and stopping powers).

The Monte Carlo engine reproduces TRIM.exe of SRIM-2013 ion by ion (same
random numbers, same Visual Basic arithmetic, same output files), and the
stopping and range tables those of SRModule.exe, digit for digit.

Simple API (energies in keV, ions and targets by name, see :mod:`lunasrim.api`)::

    import lunasrim
    rows = lunasrim.stopping_table('H', 'Si', emin=10, emax=10000)
    r = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um')
    res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=1000)
    lunasrim.material('Kapton')

Other entry points:

* :mod:`lunasrim.srim`      pysrim compatible ``Ion``, ``Layer``, ``Target``, ``TRIM``
* :func:`lunasrim.srim.run_trim`  run a :class:`lunasrim.trimin.TrimInput`
* :mod:`lunasrim.results`   readers for the output files (also of SRIM itself)
* :mod:`lunasrim.stopping`  electronic and nuclear stopping powers
* :mod:`lunasrim.sr`        stopping and range tables (SRModule.exe)
* :func:`lunasrim.straggling.energy_loss`  energy loss and straggling in a layer
* ``lunasrim`` / ``python -m lunasrim``  command line (``lunasrim --help``)

Note: ``lunasrim.trim`` is the convenience function :func:`lunasrim.api.trim`;
the TRIM engine module is imported with ``from lunasrim.trim import TrimEngine``.
"""
__version__ = '1.0.0'

from .trimin import TrimInput  # noqa: E402
from .srim import Ion, Element, Layer, Target, TRIM, run_trim  # noqa: E402
from .results import Results  # noqa: E402
from .api import Material, material, resolve_ion, stopping_table, energy_loss, trim  # noqa: E402

__all__ = ['TrimInput', 'Ion', 'Element', 'Layer', 'Target', 'TRIM', 'run_trim', 'Results',
           'Material', 'material', 'resolve_ion', 'stopping_table', 'energy_loss', 'trim',
           '__version__']
