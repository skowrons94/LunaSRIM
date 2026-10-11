"""Energy loss and energy straggling of ion beams crossing a layer (simple API).

lunasrim.energy_loss() integrates the SRIM stopping power and the
electronic energy straggling (TRIM.exe's model and Bohr's) through the
layer, without Monte Carlo.  The last case is checked against a short TRIM
run.  Output: output/energy_straggling/ (a CSV with the depth profile of
each case).  About 10 s.
"""
import csv
import os
import sys

try:
    import lunasrim
except ImportError:          # running from a checkout without installing
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    import lunasrim

OUT = os.path.abspath(os.path.join('output', 'energy_straggling'))
os.makedirs(OUT, exist_ok=True)

# ion, energy (keV), target, thickness, unit, extra arguments
CASES = [
    ('H', 400, 'Si', 1, 'um', {}),
    ('H', 100, 'C', 300, 'nm', {}),
    ('He', 2000, 'Au', 300, 'nm', {}),
    ('He', 300, 'N2', 200, 'um', {}),                  # nitrogen gas at SRIM's density
    ('C', 5000, 'Ni', 500, 'nm', {}),
    ('He', 1500, 'Mylar', 2, 'um', {'angle': 45}),       # compound of SRIM's dictionary, tilted
    ('H', 1000, 'SiO2', 100, 'ug/cm2', {'density': 2.32}),   # formula: density needed
]

print('%-4s %8s  %-6s %12s   %10s %10s %9s %9s' % ('ion', 'E (keV)', 'target', 'thickness', 'E out', 'loss',
                                                   'sigma', 'FWHM'))
for ion, energy, target, thick, unit, kw in CASES:
    r = lunasrim.energy_loss(ion, target, energy=energy, thickness=thick, unit=unit, **kw)
    print('%-4s %8g  %-6s %8g %-5s %8.2f keV %8.2f %9.3f %9.3f' % (
        ion, energy, target, thick, unit, r['energy_out'], r['energy_loss'], r['sigma'], r['fwhm']))
    name = '%s%g_%s_%g%s.csv' % (ion, energy, target, thick, unit.replace('/', '_'))
    p = r['profile']
    with open(os.path.join(OUT, name), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['depth_A', 'energy_keV', 'sigma_keV', 'sigma_bohr_keV'])
        for row in zip(p['depth'], p['energy'], p['sigma'], p['sigma_bohr']):
            w.writerow(['%.6g' % v for v in row])

# Monte Carlo check: 400 keV H through 1 um Si
res = lunasrim.trim('H', 'Si', energy=400, thickness=1, unit='um', ions=500, seed=1,
                    directory=os.path.join(OUT, 'trim_H_Si'))
e = sorted(v / 1e3 for v in res.transmit.energy)
n = len(e)
median = e[n // 2]
sigma = (e[int(0.84 * n)] - e[int(0.16 * n)]) / 2          # without the nuclear low energy tail
print('TRIM, 400 keV H in 1 um Si (%d ions): median exit energy %.1f keV, sigma %.2f keV' % (n, median, sigma))
print('CSV files in %s' % OUT)
