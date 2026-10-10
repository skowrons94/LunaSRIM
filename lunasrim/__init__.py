"""lunasrim: a Python 3 re-implementation of SRIM-2013 (TRIM and stopping powers).

The Monte Carlo engine reproduces TRIM.exe of SRIM-2013 ion by ion (same
random numbers, same Visual Basic arithmetic, same output files).

Main entry points:

* :mod:`lunasrim.srim`      pysrim compatible ``Ion``, ``Layer``, ``Target``, ``TRIM``
* :func:`lunasrim.srim.run_trim`  run a :class:`lunasrim.trimin.TrimInput`
* :mod:`lunasrim.results`   readers for the output files (also of SRIM itself)
* :mod:`lunasrim.stopping`  electronic and nuclear stopping powers
* ``python -m lunasrim TRIM.IN``  command line
"""
from .trimin import TrimInput
from .srim import Ion, Element, Layer, Target, TRIM, run_trim
from .results import Results

__version__ = '0.1.0'
__all__ = ['TrimInput', 'Ion', 'Element', 'Layer', 'Target', 'TRIM', 'run_trim', 'Results']
