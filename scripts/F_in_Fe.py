"""Implantation profile of 40 keV F in Fe."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lunasrim.srim import Ion, Layer, Target, TRIM

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

print('Average range %.1f A, straggling %.1f A' % (results.range.average_range, results.range.straggling))
