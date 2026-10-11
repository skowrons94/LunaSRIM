# LunaSRIM

lunasrim is a pure Python 3 re-implementation of SRIM-2013. It contains TRIM, the Monte Carlo transport of ions in matter, and the SR module, which tabulates stopping powers, ranges and straggling. We reverse engineered both programs from `TRIM.exe` and `SRModule.exe`, so that lunasrim reproduces them digit for digit: the same random numbers, the same Visual Basic arithmetic and the same output files. For all 46 TRIM reference cases in `tests/reference`, every output file of lunasrim is identical line by line to the file written by `TRIM.exe`, apart from the calculation date. The stopping and range tables are identical to those of `SRModule.exe` for all 24 reference cases in `tests/reference_sr`.

lunasrim needs no Windows, wine or docker. It has no required dependencies; numpy is optional and only used by the readers of the output files.

Web application: https://skowrons94.github.io/LunaSRIM/ (runs in the browser). Documentation: [docs/](https://github.com/skowrons94/LunaSRIM/blob/main/docs/index.md).

![Ion trajectories and recoil cascades of a TRIM calculation in the web application](https://raw.githubusercontent.com/skowrons94/LunaSRIM/main/docs/img/trim_trajectories.png)

## Installation

```sh
pip install lunasrim                 # the package and the lunasrim command
pip install "lunasrim[web]"          # with FastAPI and uvicorn for the local web application
pip install git+https://github.com/skowrons94/LunaSRIM    # the current version from GitHub
```

lunasrim needs Python 3.9 or later.

## Quick start

A stopping and range table of protons in silicon, from the command line:

```sh
lunasrim stopping H Si --emin 100keV --emax 1MeV
```

```text
Energy        dE/dx Elec  dE/dx Nucl    Projected Longitudinal      Lateral
                                            Range   Straggling   Straggling
---------------------------------------------------------------------------
100 keV        4.957E-01   8.673E-04       8692 A      990.7 A       1265 A
110 keV        4.836E-01   8.057E-04       9550 A       1024 A       1326 A
...
```

The stopping powers are in MeV/(mg/cm²) by default; `--units` selects any of SRIM's eight units and `--csv` or `--srim` the output format.

The energy loss and energy straggling of 400 keV protons in 1 µm of silicon, in Python:

```python
import lunasrim

r = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um')
print(r['energy_out'], r['sigma'], r['fwhm'])    # 330.46 keV, 4.48 keV, 10.56 keV
```

A TRIM calculation of 1000 argon ions of 50 keV in 200 nm of silicon, and some of its results:

```python
res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=1000, seed=1,
                    directory='ArSi')
print(res.range.average_range, res.range.straggling)    # 586.4 A, 226.9 A
print(res.vacancy.vacancies_per_ion)
```

`lunasrim.trim` writes `TRIM.IN` into `ArSi/` and TRIM's output files into `ArSi/SRIM Outputs/`. The same calculation from the command line is `lunasrim trim Ar Si --energy 50keV --thickness 200nm --seed 1`. Existing `TRIM.IN` and `SR.IN` files run with `lunasrim TRIM.IN` and `lunasrim SR.IN`.

## What it calculates

- Stopping and range tables (`lunasrim.stopping_table`, `lunasrim stopping`). Electronic and nuclear stopping, projected range and longitudinal and lateral straggling of any ion from H to U in any element, chemical formula or compound of SRIM's dictionary (356 entries), in SRIM's eight stopping units. The tables are those of `SRModule.exe`.
- Energy loss and energy straggling in a layer (`lunasrim.energy_loss`, `lunasrim eloss`). The mean exit energy and the energy spread of a beam crossing one layer, from the SRIM stopping power and TRIM's electronic straggling model, without Monte Carlo. This tool is not part of SRIM.
- TRIM Monte Carlo (`lunasrim.trim`, `lunasrim.TRIM`, `lunasrim trim`). Ion distributions, damage with the quick Kinchin-Pease estimate or full cascades, sputtering, transmitted and backscattered ions, in targets of several layers. lunasrim writes all of TRIM's output files, including `COLLISON.txt` and `EXYZ.txt`, and reads them back with `lunasrim.Results`, also when SRIM wrote them.
- Implantation profiles (web application). Gaussian depth profiles from the projected range and straggling, for one or several implantation energies and doses.

The models and their approximations are described in [docs/physics.md](https://github.com/skowrons94/LunaSRIM/blob/main/docs/physics.md).

## Validation

`tests/reference/` holds 46 `TRIM.IN` files with the output files that `TRIM.exe` wrote for them under wine. The cases cover the three calculation types, gas and compound layers, up to 8 elements, several layers, oblique and glancing incidence, plot windows, energies from 100 eV to 200 MeV, and the collision files. `tests/reference_sr/` holds 24 `SR.IN` files with the tables of `SRModule.exe` and two inputs for which `SRModule.exe` stops with a run-time error. The test suite runs lunasrim on all of them and compares the files line by line:

```sh
python -m unittest discover tests        # about 3 minutes
```

Exact agreement needs the arithmetic of the Visual Basic runtime: x87 square roots at 64-bit precision and the runtime's `Exp`, `Log`, `^`, `Atn`, `Sin` and `Cos`. lunasrim runs these instruction sequences natively on x86-64 processors that give the Intel results. Elsewhere it uses Python's `math` module, and up to one line of `HPbLong/LATERAL.txt` (two lines of `XeSiHigh/LATERAL.txt` on Windows) then differs in the last digit. [docs/validation.md](https://github.com/skowrons94/LunaSRIM/blob/main/docs/validation.md) describes the reference cases, the comparison tools and the known limits.

## Web application

The web application has tabs for stopping and range tables, implantation profiles, energy loss, TRIM with trajectories and damage plots, the comparison of stopping curves, and the physics documentation.

- Online at https://skowrons94.github.io/LunaSRIM/. lunasrim runs in your browser with Pyodide, so no calculation is done on a server. The first visit downloads about 14 MB, which the browser caches.
- Locally: `pip install "lunasrim[web]"`, then `lunasrim web` and open http://127.0.0.1:8000/.
- With Docker: `docker build -t lunasrim .` and `docker run -p 8000:8000 lunasrim`.

| | |
|---|---|
| ![Stopping and range tables](https://raw.githubusercontent.com/skowrons94/LunaSRIM/main/docs/img/stopping.png) | ![Implantation profiles](https://raw.githubusercontent.com/skowrons94/LunaSRIM/main/docs/img/implantation.png) |
| ![Energy loss and straggling](https://raw.githubusercontent.com/skowrons94/LunaSRIM/main/docs/img/energy_loss.png) | ![TRIM damage distributions](https://raw.githubusercontent.com/skowrons94/LunaSRIM/main/docs/img/trim_damage.png) |

See [docs/web.md](https://github.com/skowrons94/LunaSRIM/blob/main/docs/web.md) for the tools, the server settings and the deployment of the static site.

## Documentation

- [User guide](https://github.com/skowrons94/LunaSRIM/blob/main/docs/user-guide.md): Python API, command line, examples, performance and limits.
- [Physics](https://github.com/skowrons94/LunaSRIM/blob/main/docs/physics.md): the models of SRIM-2013 as implemented, with their approximations.
- [Validation](https://github.com/skowrons94/LunaSRIM/blob/main/docs/validation.md): the comparison with `TRIM.exe` and `SRModule.exe`.
- [Web application](https://github.com/skowrons94/LunaSRIM/blob/main/docs/web.md): using, running and deploying it.
- [Development](https://github.com/skowrons94/LunaSRIM/blob/main/docs/development.md): repository layout, tests and release process.
- [Examples](https://github.com/skowrons94/LunaSRIM/blob/main/examples/README.md): scripts for stopping tables, energy straggling and TRIM.
- [Changelog](https://github.com/skowrons94/LunaSRIM/blob/main/CHANGELOG.md).

## Citing

When you publish results obtained with lunasrim, cite SRIM, whose physics lunasrim reproduces:

J. F. Ziegler, M. D. Ziegler and J. P. Biersack, SRIM – The stopping and range of ions in matter (2010), Nucl. Instr. Meth. B 268 (2010) 1818.

To cite lunasrim itself, use the metadata in [CITATION.cff](https://github.com/skowrons94/LunaSRIM/blob/main/CITATION.cff) (GitHub shows it under "Cite this repository").

## License

The lunasrim code is released under the MIT license ([LICENSE](https://github.com/skowrons94/LunaSRIM/blob/main/LICENSE)). The SRIM-2013 data files in `lunasrim/data` (`SCOEF03.DAT`, `SNUC03.dat`, `ATOMDATA`, `Compound.dat`) are by J. F. Ziegler and remain under the terms of the SRIM authors.

## Acknowledgements

SRIM is written by J. F. Ziegler, J. P. Biersack and M. D. Ziegler (www.srim.org). The reference outputs were produced with the SRIM-2013 installation in the docker image of [pysrim](https://github.com/costrouc/pysrim) by C. Ostrouchov, whose Python interface `lunasrim.srim` follows, so that pysrim scripts run on lunasrim after changing the import.
