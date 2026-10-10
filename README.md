# LunaSRIM

A Python 3 re-implementation of **SRIM-2013**, obtained by reverse
engineering `TRIM.exe` and `SRModule.exe` from the
[costrouc/pysrim](https://github.com/costrouc/pysrim) docker image. It
covers TRIM (Monte Carlo transport of ions in matter) and the SR module
(tables of stopping powers, ranges and straggling).

TRIM is reproduced **ion by ion**: same random number sequence, same
physics, same Visual Basic arithmetic (including Single precision
variables) and the same output files. For the reference cases in
`tests/reference`, the files written by lunasrim are identical line by line
to those of the original `TRIM.exe`. The only exceptions are a few last
digits caused by math-library rounding; see [Validation](#validation).
The stopping and range tables are identical to those of `SRModule.exe`.

No Windows, wine or docker is needed to run calculations; the package is
pure Python (numpy is optional and only used by the result readers).

A **web application** runs on top of the package: stopping and range
tables, implantation profiles, TRIM Monte Carlo with trajectories and
damage plots, and comparisons of stopping curves. See
[Web application](#web-application).

![TRIM trajectories in the web application](docs/img/trim_trajectories.png)

## Installation

```sh
pip install .            # or just put the repository on PYTHONPATH
```

## Usage

### pysrim compatible API

```python
from lunasrim import Ion, Layer, Target, TRIM

ion = Ion('Ar', energy=50e3)                        # energy in eV
layer = Layer({'Si': {'stoich': 1, 'E_d': 15, 'lattice': 2, 'surface': 4.7}},
              density=2.321, width=2000)            # g/cm3, Angstrom
trim = TRIM(Target([layer]), ion, number_ions=1000, calculation=1,
            ranges=True, transmit=True, backscattered=True, sputtered=True)
results = trim.run('output/ArSi')    # TRIM.IN + 'SRIM Outputs/*.txt'

print(results.range.average_range, results.range.straggling)
print(results.vacancy.vacancies_per_ion)
depth, ioniz = results.ioniz.depth, results.ioniz.ions
```

The arguments follow pysrim's `srim.TRIM`. `calculation` is 1 (ion
distribution and quick Kinchin-Pease damage), 2 (full damage cascades) or
3 (monolayer collision steps / surface sputtering). Masses and surface
binding energies default to SRIM's values (most abundant isotope for the
ion, natural weight for target atoms). `TRIM.copy_output_files(src, dst)` is
available as in pysrim, so existing pysrim scripts only need
`from lunasrim import ...` and a directory for `run()`.

### Running a TRIM.IN file

```sh
python -m lunasrim TRIM.IN -o "SRIM Outputs"      # or: lunasrim TRIM.IN
python -m lunasrim TRIM.IN -n 100                  # override the number of ions
```

```python
from lunasrim import TrimInput, run_trim
engine = run_trim(TrimInput.read('TRIM.IN'), 'SRIM Outputs')
```

### Reading output files

`lunasrim.Results(directory)` reads RANGE, IONIZ, PHONON, E2RECOIL,
VACANCY, NOVAC, LATERAL, TRANSMIT, BACKSCAT, SPUTTER, RANGE_3D and TRIMOUT
files. It works the same for files written by the real SRIM.

### Stopping and range tables (SR module)

```sh
python -m lunasrim SR.IN            # writes the table file named in SR.IN
```

```python
from lunasrim import sr
sr.run('SR.IN')
rows = sr.stopping_table(1, 1.008, [(6, 1, 12.011)], density=2.253,
                         emin=10, emax=10000, units=7)
# rows: (E keV, S_elec, S_nuc, projected range A, long. straggling A, lat. straggling A)
```

The SR.IN format is the one of SRIM's "SR Module" and "Stopping and Range
Tables". With E-Min = E-Max = 0, a list of energies follows and only the
stopping powers are written. The range integration is the module's own
PRAL-type calculation. It runs on SRIM's energy grid with 50 sub-steps per
interval and accumulates the energy, which is why SRIM prints rows like
`999.999 eV`. For inputs where SRModule.exe stops with a run-time error (a
negative straggling variance), `sr.SrError` is raised.

The underlying stopping powers are also available directly:

```python
from lunasrim import stopping
stopping.electronic_stopping(1, 1.008, 6, 100.0)      # H in C at 100 keV, eV/(1E15 atoms/cm2)
stopping.nuclear_stopping(1, 1.008, 6, 12.011, 100.0)
stopping.compound_stopping(8, 15.999, [(1, 1.008, 2), (8, 15.999, 1)], 1000.0, phase=1)
```

### Example scripts

`scripts/` has the original LunaSRIM calculations, now running on
lunasrim:

| script | calculation |
|--------|-------------|
| `F_in_Fe.py` | implantation of 40 keV F in Fe |
| `H_in_C.py`  | energy loss of 100 keV protons in 1000 A of carbon |
| `He_in_N.py` | angular straggling of 835.8 keV He in 22 cm of N2 gas |

```sh
./run.sh scripts/F_in_Fe.py        # results in output/F_in_Fe/
```

## Web application

```sh
pip install ".[web]"                 # fastapi + uvicorn
tools/fetch_web_vendor.sh            # optional: local copy of plotly.js for offline use
python -m lunasrim.web               # http://127.0.0.1:8000/
```

Or run it with Docker: `docker build -t lunasrim . && docker run -p 8000:8000 lunasrim`.

The web application has four tools.

**Stopping & Range.** This is the SRIM "Stopping and Range Tables":
- pick the ion from a periodic table, and an element, a chemical formula or
  one of SRIM's 356 compounds (with densities) as target;
- plot electronic and nuclear stopping, projected range and straggling in any
  of SRIM's 8 stopping units;
- view the table and download it as CSV or as the SRModule.exe output file.

**Implantation.** For one or more implant energies and doses it gives:
- the projected range, longitudinal and lateral straggling, and the Gaussian
  depth profiles with their sum;
- peak concentration and atomic fraction;
- R<sub>p</sub> ± ΔR<sub>p</sub> versus energy;
- a box-profile helper that chooses energies and doses for a flat
  concentration between two depths.

**TRIM Monte Carlo.** You define multi-layer targets (elements or compounds,
gas or solid, with displacement, lattice and surface energies) and the
calculation type, then follow the run live. The results are:
- ion trajectories with recoil cascades, in three projections;
- the ion and recoil range distributions with a Gaussian fit;
- vacancies per element, replacement collisions, ionization, phonons and
  energy given to recoils;
- lateral spread and final ion positions;
- energy and angle spectra of transmitted, backscattered and sputtered atoms;
- all of TRIM's output files, also as a zip.

A depth window (SRIM's plot window) sets the histogram resolution.

**Compare.** Overlay stopping, range or straggling curves of different
ion/target pairs, versus energy or energy per nucleon.

The server limits the number of ions per TRIM calculation
(`LUNASRIM_MAX_IONS`, default 100000) and the number of calculations that run
at once (`LUNASRIM_MAX_JOBS`, default 2).

| | |
|---|---|
| ![Stopping](docs/img/stopping.png) | ![Implantation](docs/img/implantation.png) |
| ![Damage](docs/img/trim_damage.png) | ![Compare](docs/img/compare_dark.png) |

## Validation

`tests/reference/` holds `TRIM.IN` files and the output files the original
`TRIM.exe` wrote for them. The cases cover:

- Kinchin-Pease and full cascades, and sputtering mode;
- gas and compound layers;
- several layers, including non-integer widths;
- oblique incidence and a plot window.

`tests/reference_sr/` does the same for `SRModule.exe` and 13 `SR.IN`
files: H to U ions, solids and gases, all 8 stopping units, eV to GeV
energies, and Bragg corrections, in both table and list modes. All of them
are identical to the original output.

```sh
python -m unittest discover tests        # or: pytest
```

The tools to produce such references are:

```sh
tools/extract_srim.sh srim               # copy SRIM-2013 out of the docker image
tools/run_real_trim.sh some/dir          # run TRIM.exe (wine, docker) on some/dir/TRIM.IN
tools/validate.py some/dir               # run lunasrim on it and diff every output file
tools/run_real_sr.sh some/dir            # run SRModule.exe on some/dir/SR.IN
tools/validate_sr.py some/dir            # compare lunasrim's table with it
```

Known differences, in the last printed digit:

- In `HeAuC`, two of the 298 transmitted ions have a direction cosine that
  differs in the 7th digit.
- In `UFeKP`, one ion's final lateral position differs in the 5th digit.

The trajectory of that UFeKP ion was compared with TRIM.exe running in an
x86 emulator, collision by collision (3599 collisions), and the two agree
exactly. The emulator replaces the VB runtime's math functions with Python
ones. So these rare differences come from the last bits of the math library
(the Windows runtime under wine versus the C library behind Python), not
from the algorithm.

## What is implemented

- Target setup: layers, compounds, gas and solid phases, Bragg corrections,
  and SRIM's per-element electronic stopping tables.
- Ion transport: free flight path, the MAGIC scattering integral, and
  electronic energy loss with straggling.
- Recoils: Kinchin-Pease damage, full cascades, and monolayer steps with
  sputtering.
- Output files: `TDATA`, `RANGE`, `LATERAL`, `IONIZ`, `PHONON`, `E2RECOIL`,
  `VACANCY` and `NOVAC`.
- Per-ion files: `TRIMOUT`, `TRANSMIT`, `BACKSCAT`, `SPUTTER` and `RANGE_3D`.
- SR module: stopping and range tables, and stopping for energy lists.

Not implemented:

- `COLLISON.txt` and `EXYZ.txt`;
- TRIM's special modes 4-8 (neutron and ion-only calculations; partly
  present but untested);
- auto-save and restore of a calculation.

## How it was done

The binaries are Visual Basic 5 native code. They were decompiled with
Ghidra and run in a small Unicorn-based x86 emulator of the VB runtime. The
emulator traced the program's variables and random numbers, and each
Python routine was compared against it value by value. Findings that matter
for exact agreement:

- **Random numbers.** The RNG is L'Ecuyer's combined generator without
  shuffling, using VB's `CLng` (round half to even). The "Reminders" field
  of TRIM.IN is the second seed. The RNG state is restored after a recoil
  cascade (calculations 1-3), so ion trajectories do not depend on recoils.
- **Variant typing.** TRIM's variables are VB Variants, so arithmetic follows
  the OLE type rules: Integer / Integer gives a Single, Single × Integer
  gives a Single, and so on.
  - The depth bin width is therefore a Single for integer layer widths, and
    the tallies accumulate in single precision.
  - The ion mass and the recoil arrays are Singles, and the stoichiometry
    sums are done in single precision.
- **Number formatting.** VB's `Str$`/`CStr` use 15 significant digits (7 for
  Singles), while `Format$` rounds a 17-digit conversion. Some output
  columns keep a format string or a fixed-length string field from the file
  written before them.

## Repository structure

- `lunasrim/`: the package.
  - `trim.py`: the Monte Carlo engine.
  - `trim_setup.py`: target constants.
  - `stopping.py`: stopping powers.
  - `sr.py`: SR module (stopping and range tables).
  - `output.py`: output files.
  - `srim.py`: API.
  - `results.py`: readers.
  - `trimin.py`: TRIM.IN.
  - `rng.py`: random numbers.
  - `vbfmt.py`: VB number formatting.
  - `tally.py`: distributions.
  - `compounds.py`: SRIM's compound dictionary.
  - `web/`: the web application (FastAPI server and browser front end).
  - `data/`: SRIM-2013 data files (`SCOEF03.DAT`, `SNUC03.dat`, `ATOMDATA`, `Compound.dat`; (c) J. F. Ziegler).
- `scripts/`: example calculations.
- `tests/`: regression tests against TRIM.exe and SRModule.exe.
- `tools/`: helpers to extract and run the original SRIM through docker and to validate lunasrim against it.

SRIM is written by J. F. Ziegler, J. P. Biersack and M. D. Ziegler
(www.srim.org). When publishing results, cite SRIM.
