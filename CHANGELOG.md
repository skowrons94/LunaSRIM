# Changelog

All notable changes to this project are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-11

First public release. LunaSRIM is a pure Python 3 re-implementation of
SRIM-2013, reverse engineered from `TRIM.exe` and `SRModule.exe` so that it
reproduces their output files digit for digit, without Windows, wine or
docker.

### Added

- **TRIM Monte Carlo** (`TRIM.exe`), reproduced ion by ion: the same random
  number generator (L'Ecuyer, VB `CLng`), the same physics (MAGIC scattering
  integral, electronic stopping and straggling, Kinchin-Pease damage, full
  damage cascades, monolayer collision steps with sputtering), and Visual
  Basic semantics (Variant typing, Single precision tallies, `Str$`/`Format$`
  number formatting).
- **All TRIM output files**: `TDATA`, `RANGE`, `LATERAL`, `IONIZ`, `PHONON`,
  `E2RECOIL`, `VACANCY`, `NOVAC`, `TRIMOUT`, `TRANSMIT`, `BACKSCAT`,
  `SPUTTER`, `RANGE_3D`, and the collision files `COLLISON.txt` (both detail
  levels) and `EXYZ.txt`, with the same file set as TRIM. All 46 reference
  cases in `tests/reference` are identical line by line to `TRIM.exe`.
- **VB runtime x87 math** (`lunasrim.vbmath`): the `Exp`, `Log`, `^`, `Atn`,
  `Sin`, `Cos` and `Sqr` instruction sequences of msvbvm50.dll, run natively
  on x86-64 where they reproduce the Intel values, with an exact pure Python
  fallback elsewhere (only one LATERAL.txt value of one reference case then
  differs in its last digits).
- **SR module** (`SRModule.exe`, `lunasrim.sr`): stopping and range tables
  (electronic and nuclear stopping, projected range, longitudinal and lateral
  straggling) in all 8 SRIM stopping units, and stopping for energy lists;
  identical to the original for all 24 reference `SR.IN` cases.
- **Energy loss and energy straggling** through a layer
  (`lunasrim.straggling.energy_loss`): mean energy from the SRIM stopping,
  straggling from TRIM.exe's model and from Bohr's formula, thickness in
  length, mass or areal density units.
- **Simple Python API** (`lunasrim.stopping_table`, `lunasrim.energy_loss`,
  `lunasrim.trim`, `lunasrim.material`): ions and targets by name (elements,
  chemical formulas, the 356 compounds of SRIM's dictionary), SRIM's default
  densities, gas phases and binding energies, energies in keV.
- **pysrim compatible API** (`Ion`, `Layer`, `Target`, `TRIM`) and readers for
  the output files (`lunasrim.Results`), which also read files written by
  SRIM itself.
- **Command line interface** `lunasrim`: `run` (TRIM.IN and SR.IN files, also
  as `lunasrim FILE`), `stopping`, `eloss`, `trim` and `web` subcommands.
- **Web application** (`lunasrim web`, FastAPI): stopping and range tables,
  implantation profiles, energy loss and straggling, TRIM with trajectories,
  damage and spectra plots, comparison of stopping curves, and the output
  files as a zip.
- **Static website** on GitHub Pages (https://skowrons94.github.io/LunaSRIM/):
  the same application running in the browser with Pyodide.
- Example scripts (`examples/`), regression tests against `TRIM.exe` and
  `SRModule.exe` (`tests/`), validation tools (`tools/`), Docker image for the
  web application, and CI and release workflows.

[1.0.0]: https://github.com/skowrons94/LunaSRIM/releases/tag/v1.0.0
