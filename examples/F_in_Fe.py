"""Implantation profile of 40 keV F in Fe (TRIM, pysrim compatible API).

The pysrim classes Ion, Layer, Target and TRIM work as in pysrim, but
TRIM.run() executes lunasrim's Python TRIM instead of TRIM.exe.
Output: output/F_in_Fe/ (TRIM.IN and "SRIM Outputs/*.txt").  About 10 s.
"""
import os
import sys

try:
    import lunasrim  # noqa: F401
except ImportError:          # running from a checkout without installing
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lunasrim import Ion, Layer, Target, TRIM

# directory to store the output (does not have to exist)
DATA_DIRECTORY = os.path.abspath(os.path.join('output', 'F_in_Fe'))

# Construct an ion (energy in eV)
ion = Ion('F', energy=40 * 1e3)

# Construct a layer (width in Angstrom)
layer = Layer({'Fe': {'stoich': 1, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}}, density=7.87, width=5e4)

# Construct a target
target = Target([layer])

# Initialize a TRIM calculation
trim = TRIM(target, ion, number_ions=1000, calculation=1, ranges=True)

# Run it with the Python engine (no SRIM.exe / wine needed)
results = trim.run(DATA_DIRECTORY)

r = results.range
print('Average range %.1f A, straggling %.1f A' % (r.average_range, r.straggling))
print('Vacancies per ion: %.1f' % results.vacancy.vacancies_per_ion)
print('Output files in %s' % results.directory)
