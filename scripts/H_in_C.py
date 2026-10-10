"""Energy loss of 100 keV protons in 1000 A of carbon."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lunasrim.srim import Ion, Layer, Target, TRIM

DATA_DIRECTORY = os.path.abspath(os.path.join('output', 'H_in_C'))

ion = Ion('H', energy=0.1 * 1e6)
layer = Layer({'C': {'stoich': 1, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}}, density=2.2, width=1000)
target = Target([layer])

# full damage cascades; COLLISON.txt is not produced by lunasrim
trim = TRIM(target, ion, number_ions=1000, calculation=2, ranges=True, transmit=True)
results = trim.run(DATA_DIRECTORY)

t = results.transmit
mean = sum(t.energy) / len(t)
print('%d ions transmitted, mean energy %.2f keV (energy loss %.2f keV)' % (len(t), mean / 1e3, 100 - mean / 1e3))
