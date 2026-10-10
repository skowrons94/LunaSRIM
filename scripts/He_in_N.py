"""Angular straggling of 835.8 keV He in 22 cm of N2 gas."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lunasrim.srim import Ion, Layer, Target, TRIM

DATA_DIRECTORY = os.path.abspath(os.path.join('output', 'He_in_N'))

ion = Ion('He', energy=835800)
layer = Layer({'N': {'stoich': 2, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}}, phase=1, density=6.75e-7, width=2.2e9)
target = Target([layer])

trim = TRIM(target, ion, number_ions=10000, calculation=1, transmit=True)
results = trim.run(DATA_DIRECTORY)

t = results.transmit
angles = [math.degrees(math.acos(min(1.0, c))) for c in t.direction[0]]
rms = math.sqrt(sum(a * a for a in angles) / len(angles))
print('%d ions transmitted, rms angle to the beam axis %.4f deg' % (len(t), rms))
