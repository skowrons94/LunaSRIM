import os
import sys

from srim import Ion, Layer, Target, TRIM

# directory to srim executable (must be SRIM 2013)
SRIM_EXECUTABLE_DIRECTORY = os.path.abspath('/tmp/srim')
# directory to store srim output (does not have to exist)
DATA_DIRECTORY = os.path.abspath('/tmp/output')

# Construct an ion
ion = Ion('H', energy=0.1 * 1e6)

# Construct a layer
layer = Layer({'C': { 'stoich': 1, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}}, density=2.2, width=1000)

# Construct a target
target = Target([layer])

# Initialize a TRIM calculation
trim = TRIM(target, ion, number_ions=1000, calculation=2, ranges=True, transmit=True, collissions=True)

# Specify the directory of SRIM.exe
srim_executable_directory = '/tmp/srim'
results = trim.run(srim_executable_directory)

TRIM.copy_output_files(SRIM_EXECUTABLE_DIRECTORY, DATA_DIRECTORY)