"""Energy loss of 100 keV protons in 1000 A of carbon (TRIM, pysrim compatible API).

Full damage cascades (calculation=2); the energies of the transmitted
protons come from TRANSMIT.txt.  Output: output/H_in_C/.  About 3 s.
"""
import os
import sys

try:
    import lunasrim  # noqa: F401
except ImportError:          # running from a checkout without installing
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lunasrim import Ion, Layer, Target, TRIM

DATA_DIRECTORY = os.path.abspath(os.path.join('output', 'H_in_C'))

ion = Ion('H', energy=0.1 * 1e6)
layer = Layer({'C': {'stoich': 1, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}}, density=2.2, width=1000)
target = Target([layer])

# full damage cascades; collisions=True would also write COLLISON.txt, as TRIM.exe does
trim = TRIM(target, ion, number_ions=1000, calculation=2, ranges=True, transmit=True)
results = trim.run(DATA_DIRECTORY)

t = results.transmit
energies = sorted(e / 1e3 for e in t.energy)            # keV
mean = sum(energies) / len(energies)
median = energies[len(energies) // 2]
print('%d ions transmitted, mean energy %.2f keV (energy loss %.2f keV), median %.2f keV'
      % (len(t), mean, 100 - mean, median))
print('Output files in %s' % results.directory)
