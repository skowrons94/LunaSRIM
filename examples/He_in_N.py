"""Angular straggling of 835.8 keV He in 22 cm of N2 gas (simple API).

lunasrim.trim() runs TRIM for an ion in one layer, with the TRIM default
binding energies.  The gas here is at low pressure (6.75e-7 g/cm3).
lunasrim.energy_loss() gives the mean energy loss and the energy straggling
of the same layer analytically.  Output: output/He_in_N/.  About 20 s.
"""
import math
import os
import sys

try:
    import lunasrim
except ImportError:          # running from a checkout without installing
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    import lunasrim

DATA_DIRECTORY = os.path.abspath(os.path.join('output', 'He_in_N'))

# 'N2' is nitrogen (a gas element: gas=True by default); 22 cm = 2.2e9 A
res = lunasrim.trim('He', 'N2', energy=835.8, thickness=22, unit='cm', density=6.75e-7,
                    ions=10000, seed=12345, directory=DATA_DIRECTORY)

t = res.transmit
angles = [math.degrees(math.acos(min(1.0, c))) for c in t.direction[0]]
rms = math.sqrt(sum(a * a for a in angles) / len(angles))
mean_e = sum(t.energy) / len(t) / 1e3
print('%d ions transmitted, rms angle to the beam axis %.4f deg, mean energy %.2f keV' % (len(t), rms, mean_e))

r = lunasrim.energy_loss('He', 'N2', energy=835.8, thickness=22, unit='cm', density=6.75e-7)
print('analytic: exit energy %.2f keV, straggling sigma %.3f keV (%.4g atoms/cm2)'
      % (r['energy_out'], r['sigma'], r['areal_density']))
print('Output files in %s' % res.directory)
