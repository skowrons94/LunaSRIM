# LunaSRIM

A Python 3 re-implementation of **SRIM-2013** (TRIM Monte Carlo transport of
ions in matter, plus the SRIM stopping powers), obtained by reverse
engineering `TRIM.exe` and `SRModule.exe` from the
[costrouc/pysrim](https://github.com/costrouc/pysrim) docker image.

TRIM is reproduced **ion by ion**: same random number sequence, same
physics, same Visual Basic arithmetic (including Single precision
variables) and the same output files. For the reference cases in
`tests/reference`, the files written by lunasrim are identical line by line
to those of the original `TRIM.exe`. The only exceptions are a few last
digits caused by math-library rounding; see [Validation](#validation).

No Windows, wine or docker is needed to run calculations; the package is
pure Python (numpy is optional and only used by the result readers).

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

### Stopping powers

```python
from lunasrim import stopping
stopping.electronic_stopping(1, 1.008, 6, 100.0)      # H in C at 100 keV, eV/(1E15 atoms/cm2)
stopping.nuclear_stopping(1, 1.008, 6, 12.011, 100.0)
stopping.compound_stopping(8, 15.999, [(1, 1.008, 2), (8, 15.999, 1)], 1000.0, phase=1)
```

These are the SRIM-2013 stopping routines of the SR module (proton, helium
and lithium fits, heavy ion scaling, Bethe-Bloch at high energy, phase and
compound corrections). They agree with `SRModule.exe` to better than 1e-9.

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

## Validation

`tests/reference/` holds `TRIM.IN` files and the output files the original
`TRIM.exe` wrote for them. The cases cover:

- Kinchin-Pease and full cascades, and sputtering mode;
- gas and compound layers;
- several layers, including non-integer widths;
- oblique incidence and a plot window.

```sh
python -m unittest discover tests        # or: pytest
```

The tools to produce such references are:

```sh
tools/extract_srim.sh srim               # copy SRIM-2013 out of the docker image
tools/run_real_trim.sh some/dir          # run TRIM.exe (wine, docker) on some/dir/TRIM.IN
tools/validate.py some/dir               # run lunasrim on it and diff every output file
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

Not implemented:

- `COLLISON.txt` and `EXYZ.txt`;
- TRIM's special modes 4-8 (neutron and ion-only calculations; partly
  present but untested);
- auto-save and restore of a calculation;
- the SR module's range and straggling tables (only the stopping powers are
  available).

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
  - `output.py`: output files.
  - `srim.py`: API.
  - `results.py`: readers.
  - `trimin.py`: TRIM.IN.
  - `rng.py`: random numbers.
  - `vbfmt.py`: VB number formatting.
  - `tally.py`: distributions.
  - `data/`: SRIM-2013 data files (`SCOEF03.DAT`, `SNUC03.dat`, `ATOMDATA`; (c) J. F. Ziegler).
- `scripts/`: example calculations.
- `tests/`: regression tests against TRIM.exe.
- `tools/`: helpers to extract and run the original SRIM through docker and to validate lunasrim against it.

SRIM is written by J. F. Ziegler, J. P. Biersack and M. D. Ziegler
(www.srim.org). When publishing results, cite SRIM.
