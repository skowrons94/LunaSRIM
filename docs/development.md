# Development

This document describes the layout of the repository, how to run the tests and compare lunasrim with the original SRIM programs, and how a release and the website are published.

## Repository layout

- `lunasrim/`: the package.
  - `api.py`: the simple API (`stopping_table`, `energy_loss`, `trim`, `material`).
  - `__main__.py`: the command line interface.
  - `srim.py`: the pysrim-compatible `Ion`, `Layer`, `Target`, `TRIM` and `run_trim`.
  - `trim.py`: the TRIM Monte Carlo engine (`TrimEngine`).
  - `trim_setup.py`: the target constants and stopping tables that TRIM computes before the first ion.
  - `stopping.py`: electronic and nuclear stopping powers.
  - `sr.py`: the SR module (stopping and range tables, `SR.IN` files).
  - `straggling.py`: energy loss and energy straggling through a layer.
  - `trimin.py`: reading and writing `TRIM.IN`.
  - `output.py`, `tally.py`: TRIM's output files and distributions.
  - `results.py`: readers of the output files.
  - `rng.py`: TRIM's random number generator.
  - `vbmath.py`: the x87 math of the Visual Basic runtime; `vbfmt.py`: Visual Basic number formatting.
  - `compounds.py`, `data_files.py`: SRIM's compound dictionary and data files.
  - `data/`: the SRIM-2013 data files `SCOEF03.DAT`, `SNUC03.dat`, `ATOMDATA` and `Compound.dat`.
  - `web/`: the web application: `core.py` (calculations behind the API), `server.py` (FastAPI), `pyodide_api.py` (entry point in the browser), `static/` (front end, and `physics.md`, a copy of `docs/physics.md`).
- `tests/`: unit tests and the reference outputs of `TRIM.exe` (`reference/`) and `SRModule.exe` (`reference_sr/`).
- `tools/`: scripts to run the original SRIM under docker and wine, to compare with it, to build the static site and to fetch the JavaScript libraries. `vbmath.S` is the assembly source of the machine code in `vbmath.py`.
- `examples/`: example scripts.
- `docs/`: this documentation and the screenshots in `docs/img/`.
- `.github/workflows/`: `tests.yml` (tests), `pages.yml` (website) and `release.yml` (release).

## Tests

```sh
pip install -e ".[web]" httpx           # httpx for the tests of the FastAPI server
python -m unittest discover tests       # or: pytest
```

The tests compare lunasrim with the 46 `TRIM.exe` and 24 `SRModule.exe` reference cases (validation.md) and check the simple API, the command line, the web API in both forms, the analytic straggling and the x87 math. Without fastapi the tests of the server are skipped. The whole suite takes about 3 minutes on one core; `tests/test_trim.py` takes most of it. A single case runs with, for example, `python -m unittest tests.test_trim.TestAgainstTrim.test_HeArGasW`.

On a processor where `lunasrim.vbmath.NATIVE` is false (validation.md), `tests/test_trim.py` allows the one known difference in `HPbLong/LATERAL.txt`. The workflow `.github/workflows/tests.yml` runs the suite on Linux with Python 3.9, 3.10, 3.12 and 3.13 and on Windows and macOS with Python 3.12, and runs three commands of the command line.

## Comparing with the original SRIM

The original programs run under wine in the docker image `costrouc/pysrim`. A new reference case is made in four steps:

```sh
mkdir tests/reference/NewCase
cp TRIM.IN tests/reference/NewCase/              # the input to test
tools/run_real_trim.sh tests/reference/NewCase   # TRIM.exe writes its output files there
tools/validate.py tests/reference/NewCase        # lunasrim on the same input, file by file
```

Remove `TRIM.IN.used` and any file that is not an output file before committing the case; `tests/test_trim.py` picks up every directory of `tests/reference/` and requires the same set of `.txt` files. For the SR module, `tools/run_real_sr.sh` and `tools/validate_sr.py` work in the same way with `SR.IN` and `tests/reference_sr/`. `tools/extract_srim.sh` copies the SRIM-2013 installation out of the image, for example to look at the binaries or to refresh `lunasrim/data`.

A change of the physics or of the arithmetic must keep all reference cases identical. When a case differs, the first differing line and the trajectory of that ion usually point to the routine; a one-ion TRIM.IN with the random number state of that ion as its two seeds reproduces it in isolation.

## Documentation

The documentation is in `docs/`. `docs/physics.md` is also shown on the Physics tab of the web application, from the packaged copy `lunasrim/web/static/physics.md`. After editing it, run

```sh
python tools/sync_docs.py            # copy docs/physics.md into the package
python tools/sync_docs.py --check    # exit status 1 if the copy is out of date
```

`tests/test_docs_page.py` fails while the two files differ. The formulas of `physics.md` are rendered with KaTeX on the website and by GitHub, so they use plain LaTeX without custom macros.

## Releases

A release is published by pushing a version tag:

1. Set `__version__` in `lunasrim/__init__.py` to the new version, for example `1.0.1`. `pyproject.toml` reads the version from there.
2. Add a section `## [1.0.1] - YYYY-MM-DD` to `CHANGELOG.md`, and update `version` and `date-released` in `CITATION.cff`.
3. Commit, tag and push:

   ```sh
   git commit -am "Release 1.0.1"
   git tag v1.0.1
   git push origin main v1.0.1
   ```

The workflow `.github/workflows/release.yml` then:

- checks that the tag equals `v` followed by `lunasrim.__version__`, fetches the JavaScript libraries for the offline web application, builds the source distribution and the wheel, and tests the wheel in a fresh virtual environment;
- creates the GitHub release for the tag with the distribution files, and the section of `CHANGELOG.md` for this version as release notes;
- uploads the distribution to PyPI by trusted publishing, without an API token.

Trusted publishing needs a one-time setup on PyPI. In the account on pypi.org, under Publishing (for a new project, a "pending publisher"), add a GitHub publisher with the owner `skowrons94`, the repository `LunaSRIM`, the workflow `release.yml` and the environment `pypi`. In the GitHub repository settings, create the environment `pypi` under Environments; protection rules such as required reviewers are optional.

## Website

The workflow `.github/workflows/pages.yml` builds the static site with `tools/build_static.py` and deploys it to GitHub Pages on every push to `main` (web.md). It needs the one-time setting Pages, Source: "GitHub Actions" in the repository settings.
