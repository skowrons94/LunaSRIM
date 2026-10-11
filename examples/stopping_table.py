"""Stopping and range tables (simple API, identical to SRIM's SRModule.exe).

Writes the table of H in Si as CSV and in the SRModule.exe file format,
and prints a few rows of other ion / target pairs in other units.
Output: output/stopping_table/.  About 2 s.
"""
import csv
import os
import sys

try:
    import lunasrim
except ImportError:          # running from a checkout without installing
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    import lunasrim
from lunasrim import api

OUT = os.path.abspath(os.path.join('output', 'stopping_table'))
os.makedirs(OUT, exist_ok=True)

# H in Si, 10 keV - 10 MeV, MeV/(mg/cm2)
rows = lunasrim.stopping_table('H', 'Si', emin=10, emax=10000)
with open(os.path.join(OUT, 'H_in_Si.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['energy_keV', 'S_elec_MeV/(mg/cm2)', 'S_nuc_MeV/(mg/cm2)', 'range_A', 'long_straggling_A',
                'lat_straggling_A'])
    w.writerows(rows)
# the same table in the format of SRModule.exe
api.stopping_report('H', 'Si', emin=10, emax=10000).write(os.path.join(OUT, 'Hydrogen in Si.txt'))

print('H in Si:  E (keV)   Se, Sn MeV/(mg/cm2)    range (um)')
for e, se, sn, rp, lon, lat in rows[::13]:
    print('       %9.4g   %9.4g  %9.3g   %9.4g' % (e, se, sn, rp / 1e4))

# maximum of the electronic stopping (Bragg peak) of He in Kapton, keV/um
m = lunasrim.material('Kapton')
rows = lunasrim.stopping_table('He', m, emin=10, emax=10000, units='keV/um')
e, se = max(((r[0], r[1]) for r in rows), key=lambda x: x[1])
print('He in %s (%.3g g/cm3): maximum stopping %.1f keV/um at %.0f keV' % (m.name, m.density, se, e))

# a formula needs its density; stopping in eV/(1E15 atoms/cm2)
rows = lunasrim.stopping_table('He', 'SiO2', density=2.32, emin=2000, emax=2000, units='eV/(1E15 atoms/cm2)')
print('2 MeV He in SiO2: %.2f eV/(1E15 atoms/cm2)' % rows[-1][1])
print('Files in %s' % OUT)
