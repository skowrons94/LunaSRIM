# Examples

Each script runs in less than a minute and writes its files under
`./output/<name>/` (relative to the current directory):

```sh
pip install lunasrim            # or: pip install .  in the repository
python examples/stopping_table.py
```

| script | calculation | API |
|--------|-------------|-----|
| `stopping_table.py` | stopping and range tables (H in Si, He in Kapton, He in SiO2), CSV and SRModule format | simple |
| `energy_straggling.py` | energy loss and energy straggling of several beams through layers, with a TRIM check | simple |
| `He_in_N.py` | angular straggling of 835.8 keV He in 22 cm of N2 gas (TRIM) | simple |
| `F_in_Fe.py` | implantation profile of 40 keV F in Fe (TRIM) | pysrim compatible |
| `H_in_C.py` | energy loss of 100 keV protons in 1000 A of carbon, full cascades (TRIM) | pysrim compatible |

The same calculations are available from the command line, e.g.
`lunasrim stopping H Si`, `lunasrim eloss H Si --energy 400keV --thickness 1um`
and `lunasrim trim F Fe --energy 40keV --thickness 5um`.
