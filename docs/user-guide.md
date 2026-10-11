# User guide

This guide describes how to install lunasrim and how to use it from Python and from the command line. The physical models behind the numbers are in [physics.md](physics.md), and the web application is described in [web.md](web.md).

## Installation

```sh
pip install lunasrim                 # package and the lunasrim command
pip install "lunasrim[web]"          # also FastAPI and uvicorn, for the local web application
pip install "lunasrim[numpy]"        # also numpy, for the readers of the output files
pip install git+https://github.com/skowrons94/LunaSRIM     # current version from GitHub
```

From a checkout of the repository, `pip install .` installs the package, and `pip install -e .` installs it in editable mode. lunasrim needs Python 3.9 or later and has no required dependencies. When numpy is installed, the readers of `lunasrim.results` return numpy arrays; otherwise they return lists. The SRIM-2013 data files (`SCOEF03.DAT`, `SNUC03.dat`, `ATOMDATA`, `Compound.dat`) are part of the package.

`lunasrim --version` prints the installed version.

## Python API

The functions `lunasrim.stopping_table`, `lunasrim.energy_loss`, `lunasrim.trim` and `lunasrim.material` take ions and targets by name, energies in keV and thicknesses with a unit. They resolve names and defaults and then call the low-level functions described further below, so the results are the same.

### Ions

An ion is an element symbol (`'He'`), an element name (`'Helium'`) or an atomic number (`2`), from H to U. Its mass defaults to the mass of the most abundant isotope, as in SRIM (4.003 u for He); the argument `mass=` sets another value in u. `lunasrim.resolve_ion('He')` returns the atomic number and mass that the functions will use.

### Targets and materials

A target is resolved by `lunasrim.material`, which returns a `lunasrim.api.Material`. The text is tried in this order:

1. an element symbol (`'Si'`), element name (`'silicon'`) or atomic number;
2. the name of a compound of SRIM's dictionary (`Compound.dat`, 356 entries), with or without its reference in parentheses: `'Mylar'`, `'Water_Liquid'`, `'Kapton Polyimide Film (ICRU-179)'`;
3. a chemical formula: `'SiO2'`, `'Al2O3'`, `'C8H8'`; a formula of a single element, such as `'N2'`, is that element;
4. the beginning of a compound name or description, or words of it, if they identify one entry: `'Kapton'`, `'quartz'`.

A dictionary such as `{'Si': 1, 'O': 2}` (symbol or atomic number to stoichiometry) is also accepted. An ambiguous text gives a `ValueError` that lists the candidates:

```python
>>> lunasrim.material('Kapton')
<Material 'Kapton Polyimide Film (ICRU-179)': H0.25641 C0.564103 N0.0512821 O0.128205, 1.42 g/cm3, solid>
>>> lunasrim.material('SiO2', density=2.32)
<Material 'SiO2': Si1 O2, 2.32 g/cm3, solid>
>>> lunasrim.material('N2')
<Material 'N2': N2, 0.001251 g/cm3, gas>
>>> lunasrim.material('Water')
ValueError: 'Water' matches several compounds of the SRIM dictionary: 'Water_Liquid (ICRU-276)' (1 g/cm3); ...
```

The defaults are:

- `density` (g/cm³): for an element, the solid density of SRIM's tables, or SRIM's gas density for H, He, N, O, F, Ne, Cl, Ar, Br, Kr, I, Xe, At and Rn; for a compound, its dictionary density. A formula with several elements and a dictionary need `density=`.
- `gas`: true for these gaseous elements and for gas compounds of the dictionary, false otherwise. It selects SRIM's gas-phase stopping (physics.md, Section 2.4).
- `correction`: the compound (Bragg) correction factor, 1.
- Atomic masses are the natural atomic weights.

The attributes of a `Material` are `name`, `elements` (a list of `(Z, stoichiometry, mass)`), `density`, `gas`, `correction`, `mean_mass` and `atom_density` (atoms/cm³). A `Material` can be passed wherever a target is expected.

### Stopping and range tables

```python
rows = lunasrim.stopping_table('He', 'Au', emin=1000, emax=2000, units='keV/um')
energy, s_elec, s_nuc, rp, long_strag, lat_strag = rows[0]
# 1000.0 keV, 731.8 keV/um, 1.554 keV/um, 15088 A, 2990 A, 3426 A
```

`stopping_table(ion, target, emin=10, emax=10000, units='MeV/(mg/cm2)', density=None, gas=None, mass=None, correction=None)` returns one tuple per energy of SRIM's energy grid between `emin` and `emax` (keV, or a string such as `'2MeV'`). Each tuple holds the energy (keV), the electronic and nuclear stopping in `units`, and the projected range, longitudinal straggling and lateral straggling in Å. The table is that of `SRModule.exe`. The units are 1 to 8 or one of the names `'eV/A'`, `'keV/um'`, `'MeV/mm'`, `'keV/(ug/cm2)'`, `'MeV/(mg/cm2)'`, `'keV/(mg/cm2)'`, `'eV/(1E15 atoms/cm2)'` and `'LSS'` (reduced units).

`lunasrim.api.stopping_report` takes the same arguments and returns the `lunasrim.sr.StoppingTable`. Its method `write(path)` writes the file in the format of `SRModule.exe`, and `lines()` returns its lines.

### Energy loss and energy straggling

```python
r = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um')
r['energy_out'], r['energy_loss'], r['sigma'], r['fwhm']    # 330.46, 69.54, 4.48, 10.56 keV
```

`energy_loss(ion, target, energy, thickness, unit='A', angle=0, density=None, gas=None, mass=None, correction=None)` integrates the mean energy and the energy variance of a beam through one layer (physics.md, Section 6). The thickness unit is one of `'A'`, `'nm'`, `'um'`, `'mm'`, `'cm'`, `'ug/cm2'`, `'mg/cm2'`, `'1e15/cm2'` and `'atoms/cm2'`, or the thickness is a string with its unit such as `'1um'`. `angle` tilts the layer: the path length is the thickness divided by the cosine of the angle (degrees). The result is a dictionary:

| key | meaning |
|---|---|
| `energy_out` | mean exit energy (keV) |
| `energy_loss`, `electronic_loss`, `nuclear_loss` | mean energy loss and its electronic and nuclear parts (keV) |
| `mean_energy` | mean energy averaged over the path (keV) |
| `sigma`, `fwhm` | energy straggling, one standard deviation and FWHM (keV), TRIM's model |
| `sigma_bohr`, `fwhm_bohr` | the same with Bohr's formula |
| `thickness` | path length (Å) |
| `areal_density`, `mass_thickness`, `atom_density` | atoms/cm², µg/cm² and atoms/cm³ |
| `stopped`, `depth_reached` | true if the mean energy falls below 1 % of the initial energy in the layer, and the depth reached (Å) |
| `mean_z`, `trim_z` | stoichiometric mean atomic number, and TRIM's integer mean |
| `profile` | lists along the path: `depth` (Å), `energy`, `sigma`, `sigma_bohr` (keV), `se`, `sn` (eV/Å) |
| `angle`, `path_factor`, `material` | the tilt, 1/cos(angle) and the `Material` used |

The spread includes the electronic straggling only. Nuclear collisions give a low-energy tail and path-length differences broaden the spectrum further; a TRIM calculation shows both.

### TRIM calculations in one layer

```python
res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=1000, seed=1, directory='ArSi')
res.range.average_range, res.range.straggling       # 586.4 A, 226.9 A
```

`trim(ion, target, energy, thickness, unit='A', ions=1000, calculation=1, angle=0, seed=None, directory=None, density=None, gas=None, mass=None, correction=None, progress=None)` runs TRIM for a beam in a single layer. `calculation` selects the calculation type:

1. ion distribution and quick Kinchin-Pease damage;
2. full damage cascades;
3. monolayer collision steps for surface sputtering.

`seed` is the random number seed; the same seed gives the same ions, and `None` draws a random one. The target atoms get the defaults of the TRIM user interface: displacement energies of 10 eV for H, 28 eV for C, N and O, 15 eV for Si and 25 eV for all other elements, a lattice binding energy of 3 eV, and the heat of sublimation of SRIM's tables as surface binding energy. `RANGE_3D.txt`, `BACKSCAT.txt`, `TRANSMIT.txt` and `SPUTTER.txt` are switched on. The function writes `TRIM.IN` into `directory` and the output files into `directory/SRIM Outputs`; with `directory=None` it makes a temporary directory and keeps it. `progress(ion_number, engine)` is called after every ion. The return value is a `lunasrim.Results`.

For several layers, other binding energies or the collision files, use the pysrim-compatible classes below.

### pysrim-compatible API

The classes `Ion`, `Layer`, `Target` and `TRIM` follow those of [pysrim](https://github.com/costrouc/pysrim). A pysrim script runs on lunasrim after changing the import to `from lunasrim import Ion, Layer, Target, TRIM` and giving `run()` a directory. `TRIM.run()` executes lunasrim instead of `TRIM.exe`.

```python
from lunasrim import Ion, Layer, Target, TRIM

ion = Ion('He', energy=5e6)                                   # energy in eV
mylar = Layer.from_formula('C10H8O4', density=1.397, width=20000)            # width in A
air = Layer({'N': 0.784, 'O': 0.211, 'Ar': 0.005}, density=0.00120484, width=5e6, phase=1)
si = Layer({'Si': {'stoich': 1, 'E_d': 15, 'lattice': 2, 'surface': 4.7}}, density=2.3212, width=1e6)
trim = TRIM(Target([mylar, air, si]), ion, number_ions=1000, calculation=1, random_seed=7,
            ranges=True, transmit=True, backscattered=True)
results = trim.run('HeMylarAirSi')          # TRIM.IN and 'SRIM Outputs/*.txt'
```

- `Ion(identifier, energy, mass=None)`: element symbol, name or atomic number; `energy` in eV.
- `Layer(elements, density, width, phase=0, name=None, bragg_correction=1.0)`: `elements` maps an element to its stoichiometry, or to a dictionary with `stoich` and optionally `E_d` (displacement energy), `lattice` and `surface` (binding energies), all in eV. The defaults are 25 eV, 3 eV and the heat of sublimation. `density` is in g/cm³, `width` in Å, and `phase` is 0 for a solid and 1 for a gas. `Layer.from_formula(formula, density, width, phase=0)` builds a layer from a formula.
- `Target(layers)`: the layers in the order in which the ion meets them.
- `TRIM(target, ion, number_ions=1000, calculation=1, angle_ions=0, random_seed=None, reminders=0, plot_xmin=0, plot_xmax=0, ranges=False, backscattered=False, transmit=False, sputtered=False, collisions=False, exyz=0, bragg_correction=1.0, description=...)`: `calculation` as above; `random_seed` and `reminders` are the two seeds of the random number generator; `plot_xmin` and `plot_xmax` (Å) set the depth window of the distributions (0 and 0 for the whole target); the flags switch on `RANGE_3D.txt`, `BACKSCAT.txt`, `TRANSMIT.txt` and `SPUTTER.txt`; `collisions` is 0, 1 (ion collisions) or 2 (ion collisions and the listing of every cascade) for `COLLISON.txt`; `exyz` is the energy interval (eV) of `EXYZ.txt`, 0 for none.
- `TRIM.run(directory, progress=None)` writes `TRIM.IN` into `directory`, runs the calculation and returns `Results(directory + '/SRIM Outputs')`. `TRIM.copy_output_files(src, dst)` copies the files as in pysrim.

### Reading results

`lunasrim.Results(directory)` reads all TRIM output files found in a `SRIM Outputs` directory, also when SRIM itself wrote them. Missing files give `None`. The attribute names follow pysrim:

| attribute | file | content |
|---|---|---|
| `range` | `RANGE.txt` | `depth`, `ions` (atoms/cm³ per atoms/cm²), `elements` (recoils), `average_range`, `straggling`, `lateral_range`, `lateral_straggling`, `radial_range`, `radial_straggling` |
| `ioniz` | `IONIZ.txt` | `depth`, `ions`, `recoils` (eV/(Å ion)) |
| `phonons` | `PHONON.txt` | `depth`, `ions`, `recoils` (eV/(Å ion)) |
| `e2recoil` | `E2RECOIL.txt` | `depth`, `ions`, `absorbed` |
| `vacancy` | `VACANCY.txt` | `depth`, `knock_ons`, `vacancies`, `vacancies_per_ion`, `displacements_per_ion`, `replacements_per_ion` |
| `novac` | `NOVAC.txt` | `depth`, `number` (replacement collisions) |
| `lateral` | `LATERAL.txt` | `depth`, `lateral_range`, `lateral_straggling`, `radial_range`, `radial_straggling` |
| `transmit`, `backscat`, `sputter`, `trimout` | `TRANSMIT.txt`, ... | `ion_number`, `z`, `energy` (eV), `position` (x, y, z in Å), `direction` (cosines); `len()` is the number of atoms |
| `range3d` | `RANGE_3D.txt` | `ion_number`, `positions` (x, y, z in Å) of the stopped ions |
| `exyz` | `EXYZ.txt` | `ion_number`, `energy` (keV), `positions`, `se` (eV/Å), `recoil_energy` (eV); `ion(n)` gives the rows of ion n |
| `collisions` | `COLLISON.txt` | `rows`: one dictionary per collision with energy, position, atom, recoil energy and the damage it caused |

The distribution files also have `num_ions`, the raw `lines` and the columns in `data`.

### Running SRIM input files

A `TRIM.IN` file written by SRIM or by hand runs with

```python
from lunasrim import TrimInput, run_trim
engine = run_trim(TrimInput.read('TRIM.IN'), 'SRIM Outputs')     # optional: number_ions=100
```

`run_trim` writes the same files as `TRIM.exe` and returns the engine; `engine.tal.n_back` and `engine.tal.n_trans` are the numbers of backscattered and transmitted ions. An `SR.IN` file of SRIM's "Stopping and Range Tables" runs with

```python
from lunasrim import sr
path = sr.run('SR.IN')                  # writes the table file named in SR.IN next to SR.IN
table = sr.StoppingTable.from_file('SR.IN')
```

With E-Min = E-Max = 0 the SR.IN file holds a list of energies and only the stopping powers are written. `sr.run` raises `sr.SrError` where `SRModule.exe` stops with a run-time error.

### Low-level functions

The simple API calls these functions, which take atomic numbers and masses directly:

```python
from lunasrim import sr, stopping, straggling

# stopping and range table: ion Z and mass, [(Z, stoichiometry, mass)], density g/cm3, keV, units 1-8
rows = sr.stopping_table(1, 1.008, [(6, 1, 12.011)], 2.253, 10, 10000, phase=0, correction=1, units=7)

# energy loss and straggling, thickness in unit ('A', 'nm', 'um', 'mm', 'ug/cm2', 'mg/cm2', '1e15/cm2', 'atoms/cm2')
r = straggling.energy_loss(1, 1.008, 400, [(14, 1, 28.086)], 2.3212, 1, unit='um')

# stopping cross sections in eV/(1E15 atoms/cm2), energies in keV
stopping.electronic_stopping(1, 1.008, 6, 100.0)                  # H in C at 100 keV
stopping.nuclear_stopping(1, 1.008, 6, 12.011, 100.0)
stopping.compound_stopping(8, 15.999, [(1, 1.008, 2), (8, 15.999, 1)], 1000.0, phase=1)   # (S_e, S_n)
```

## Command line

The command `lunasrim` (or `python -m lunasrim`) has five subcommands. Energies are in keV, or a number with eV, keV, MeV or GeV (`400`, `400keV`, `2MeV`). Thicknesses are a number with A, nm, um, mm, cm, ug/cm2, mg/cm2, 1e15/cm2 or atoms/cm2 (`1um`, `200nm`, `2e19atoms/cm2`). `lunasrim COMMAND -h` lists the options of a command.

### Input files: `lunasrim run`

```sh
lunasrim TRIM.IN                     # TRIM calculation, output files in "SRIM Outputs"
lunasrim TRIM.IN -n 100 -o out       # 100 ions instead of the number in TRIM.IN, output in out/
lunasrim SR.IN                       # stopping table, written next to SR.IN as named in it
lunasrim run case1/TRIM.IN case2/SR.IN    # several files; TRIM output next to each file
```

`lunasrim FILE` is short for `lunasrim run FILE`; the type of file is recognised from its first line. Options: `-o/--output DIR` (one file only), `-n/--number-ions N`, `-q/--quiet`.

### Stopping tables: `lunasrim stopping`

```sh
lunasrim stopping H Si                                     # 10 keV to 10 MeV, MeV/(mg/cm2)
lunasrim stopping He Kapton --emin 1MeV --emax 2MeV --units keV/um --csv
lunasrim stopping O SiO2 --density 2.32 --srim -o OSiO2.txt   # file in the format of SRModule.exe
```

Options: `--emin` (default `10keV`), `--emax` (default `10MeV`), `--units` (1 to 8 or a name, default `MeV/(mg/cm2)`), `--csv` or `--srim` (output format), `-o/--output FILE`, and the target options below.

### Energy loss: `lunasrim eloss`

```sh
lunasrim eloss H Si --energy 400keV --thickness 1um
```

```text
Ion:    H (Z=1, 1.008 amu), 400 keV
Target: Si [Si1], 2.3212 g/cm3, solid
Layer:  1um = 1 um
        4.977e+18 atoms/cm2, 232.12 ug/cm2

exit energy           330.464 keV
energy loss            69.536 keV  (electronic 69.4652, nuclear 0.0707167)
mean energy           365.773 keV  (average over the layer)
straggling sigma       4.4826 keV  FWHM 10.556 keV  (TRIM model)
straggling sigma       4.4765 keV  FWHM 10.541 keV  (Bohr)
```

Options: `--energy` and `--thickness` (required), `--angle DEG` (tilt of the layer), `--json` (all results as JSON) with `--profile` (include the depth profile), and the target options.

### TRIM: `lunasrim trim`

```sh
lunasrim trim Ar Si --energy 50keV --thickness 200nm --ions 1000 --seed 1 --out ArSi
```

Options: `--energy` and `--thickness` (required), `--ions` (default 1000), `--calculation 1|2|3` (default 1), `--angle DEG`, `--seed`, `--out DIR` (default `lunasrim_output`), `-q/--quiet`, and the target options. The command prints the mean range, straggling, lateral range, the numbers of backscattered and transmitted ions and the vacancies per ion, and leaves `TRIM.IN` and `SRIM Outputs/` in the output directory.

### Target options

`stopping`, `eloss` and `trim` take an ion and a target as described above, and the options `--density G_CM3`, `--gas` or `--solid`, `--mass AMU` (ion mass) and `--correction` (compound correction).

### Web application: `lunasrim web`

`lunasrim web [--host 127.0.0.1] [--port 8000]` starts the web application (requires `pip install "lunasrim[web]"`). See [web.md](web.md).

## Examples

The directory `examples/` contains complete scripts. Each writes its files under `./output/<name>/`.

| script | calculation | API |
|---|---|---|
| `stopping_table.py` | stopping and range tables (H in Si, He in Kapton, He in SiO2), as CSV and in the format of `SRModule.exe` | simple |
| `energy_straggling.py` | energy loss and straggling of several beams through layers, with a TRIM check | simple |
| `He_in_N.py` | angular straggling of 835.8 keV He in 22 cm of N₂ gas (TRIM) | simple |
| `F_in_Fe.py` | implantation profile of 40 keV F in Fe (TRIM) | pysrim compatible |
| `H_in_C.py` | energy loss of 100 keV protons in 1000 Å of carbon with full cascades (TRIM) | pysrim compatible |

## Performance

lunasrim is pure Python, and a TRIM calculation runs on one processor core. We measured these times on an Intel Xeon processor at 2.80 GHz with Python 3.13:

| calculation | time |
|---|---|
| stopping table, H in Si, 10 keV to 10 MeV | 0.13 s |
| energy loss and straggling, 400 keV H in 1 µm Si | 0.03 s |
| TRIM, 1000 ions of 400 keV H through 1 µm Si (transmitted), type 1 | 1.2 s |
| TRIM, 1000 ions of 400 keV H stopped in 6 µm Si, type 1 | 2.4 s |
| TRIM, 1000 ions of 50 keV Ar in 200 nm Si, type 1 | 6.8 s |
| TRIM, 1000 ions of 50 keV Ar in 200 nm Si, type 2 (full cascades) | 360 s |

The time of a TRIM calculation grows with the number of collisions per ion. Full cascades follow every recoil and are slower than the quick Kinchin-Pease calculation by a factor that grows with the ion mass and energy; for 50 keV Ar in Si the factor is 53. The x87 routines that reproduce SRIM's arithmetic are called through ctypes. With the pure Python math (`LUNASRIM_NO_NATIVE_MATH=1`) the Ar run of the table took 6.5 s instead of 7.3 s when both were started from the command line. The web application in the browser (Pyodide) is slower than a local Python installation.

## Limits

- Ions and target elements from Z = 1 to Z = 92.
- Energies: the SR module tabulates from 1 eV. SRIM's phase correction overflows above 5.92 GeV/u for protons and 5.85 GeV/u for heavier ions, where `SRModule.exe` and `TRIM.exe` stop with a run-time error. `lunasrim.sr.run` then raises `SrError`, a TRIM calculation `VBOverflow` and `energy_loss` a `StragglingError`; `stopping_table` ends the table at the last grid energy below the limit. The web application accepts ion energies up to 10 GeV.
- Targets are amorphous. Crystal structure, channelling, temperature and changes of the target during the irradiation are not modelled (physics.md, Section 5).
- TRIM calculation types 4 to 8 (ion-only and neutron calculations) are untested, and auto-save and restart are not implemented.
- The web application limits a TRIM calculation to 100000 ions and, as a server, runs at most 2 calculations at a time by default (web.md).
- The analytic energy straggling neglects nuclear and path-length straggling (physics.md, Section 6).
