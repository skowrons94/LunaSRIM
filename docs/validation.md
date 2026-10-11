# Validation against SRIM-2013

lunasrim is meant to give the same results as SRIM-2013, not results that agree with it within some tolerance. We therefore validate it by running the original programs and lunasrim on the same input files and comparing the output files line by line. This document describes how the original programs are run, which reference cases are kept in the repository, what "identical" means in practice, and where lunasrim and SRIM can differ.

## Running the original programs

The SRIM-2013 installation of the docker image [costrouc/pysrim](https://github.com/costrouc/pysrim) contains `TRIM.exe` and the SR module `SRModule.exe` (file version SRIM-2012.01). Both are Visual Basic 5 programs. The scripts in `tools/` run them under wine in that image:

```sh
tools/extract_srim.sh srim          # copy the SRIM-2013 installation out of the image
tools/run_real_trim.sh CASE         # run TRIM.exe on CASE/TRIM.IN, copy the output files into CASE
tools/validate.py CASE              # run lunasrim on CASE/TRIM.IN and compare every output file
tools/run_real_sr.sh CASE           # run SRModule.exe on CASE/SR.IN
tools/validate_sr.py CASE           # compare lunasrim's table with the one of SRModule.exe
```

`run_real_trim.sh` starts `TRIM.exe` in batch mode (`TRIMAUTO` set to 1) under `xvfb-run` and wine, with a time limit of 3000 s that the variable `TRIM_TIMEOUT` changes. `run_real_sr.sh` does the same for `SRModule.exe` with a limit of 600 s (`SR_TIMEOUT`). The validation scripts exit with status 1 if any file differs and print the first differing lines.

## What "identical" means

Two output files are identical when all their lines are equal as byte strings (Latin-1 text with CR LF line ends), except the line with the calculation date and time. For TRIM the set of files must also be the same: lunasrim writes exactly the files that `TRIM.exe` writes for the switches of the TRIM.IN file, and no others. This includes the per-ion files `TRIMOUT.txt`, `TRANSMIT.txt`, `BACKSCAT.txt`, `SPUTTER.txt` and `RANGE_3D.txt` and the collision files `COLLISON.txt` and `EXYZ.txt`.

A TRIM output file depends on every trajectory. A single different random number, a different number of iterations in the scattering integral or a flight path cut differently at a layer boundary changes the following trajectory and, in general, a printed digit. Identical files therefore test the ion transport ion by ion and collision by collision, not only the average results. The reference cases use 5 to 500 ions each so that they run in a few minutes.

## TRIM reference cases

`tests/reference/` contains 46 cases: 23 with the quick Kinchin-Pease calculation (type 1), 16 with full damage cascades (type 2) and 7 with monolayer collision steps for sputtering (type 3). Each directory holds the `TRIM.IN` file and the output files of `TRIM.exe`; `tests/reference/README.md` lists them in a table. Grouped by what they test:

- Basic transport and damage: `ArSi` (50 keV Ar in Si), `AuSiKP` (1 MeV Au in Si), `ArSiFC` (recoil cascades), `HW` (10 keV H on W at 45°, backscattering), `HeAuC` (Au on C, transmission).
- Gas layers: `HeN` (835.8 keV He in a 22 cm N₂ layer), `NTaN` and `NTaNFC` (TaN with a N₂ gas layer at 60°), `HeArGasW` (2 MeV He through 1 mm of Ar gas into W), `HeMylarAir` (Mylar, 0.5 mm of air and Si, with missing surface energies replaced by the defaults).
- Compounds: `XeGaAs` (four elements in two compounds), `UFeKP` (50 MeV U in steel with fractional stoichiometries), `HKaptonBR` (Kapton), `OBoneFC` (bone, six elements), `CSteel8` (an 8-element steel), `WebSiO2`, `WebAirPMMA` and `WebKaptonFC` (compounds from SRIM's dictionary as the web application builds them, up to 7 elements), `HeBragg` (compound corrections 0.92 and 1.07), `SiSiO2Si` (the same element in three layers).
- Sputtering mode: `ArCuSP`, `ArGaNSP` (2 keV Ar on GaN at 45°), `ArWLow` (100 eV Ar on W), `NAlGaAsSP`, `UAgFCSP` (50 keV U on Ag and AgAu with full cascades).
- Geometry and binning: `OAlFC` (three layers with non-integer widths), `BSiPW` and `LiNiWin` (plot windows, one from 5000 to 17500.5 Å), `NeSiGlance` (1 keV Ne on Si at 85°, mostly backscattered).
- Energy extremes and mass ratios: `AuCKP` (2 MeV Au in graphite), `HPbLong` (10 MeV H through 0.4 mm of Pb), `XeSiHigh` (200 MeV Xe in Si).
- Collision files: `ColKP`, `ColKPcomp`, `ColFC1`, `ColFC1comp`, `ColFC`, `ColFC2Lay`, `ColSP`, `ColSPcomp`, `ColGas2L` (`COLLISON.txt` at both detail levels, with and without the listing of every cascade), `ExyzKP`, `ExyzFCtr`, `ExyzBack` (`EXYZ.txt` with stopped, transmitted and backscattered ions).
- Output switches: `C02backscatonly` and `C04sputteronly` (only one per-ion file switched on, which changes the content of `TRIMOUT.txt`).

All 46 cases are identical, without exceptions, when lunasrim uses its native x87 routines (see below).

## SR module reference cases

`tests/reference_sr/` contains 26 `SR.IN` files. For 24 of them the directory also holds the file `SRModule.exe` wrote, and lunasrim's file is identical to it:

- Tables of elements: `HC` (H in C, 1 keV to 10 MeV), `HeAu`, `FeFe` (E-Min above E-Max), `USi` (U in Si up to 10 GeV), `ArAr` (Ar gas).
- Tables of compounds: `CBone` (6 elements, correction 0.97), `CPE` (polyethylene up to 1 GeV), `LiLiF`, `NeGaAs`, `NeSteel` (fractional stoichiometries), `OSiCLow` (from 1.1 keV), `OSiO2` (correction 1.05), `HProp` (1,3-propanediol, correction 0.9457121), `HWaterCorr` (water, correction 0.92), `HH2Ocorr0` (correction 0, which `SRModule.exe` treats as 1).
- Gas targets: `ArAirGas` (air), `AuCO2Gas` (CO₂), `UH2gas` (H₂).
- Lists of energies (E-Min = E-Max = 0): `HList`, `HeKaptonList`, `HePMMAList`, `LiSteelList`, `OH2O` (water vapour, 13 energies from 0.5 keV to 100 MeV, correction 0.9) and `XeList`.

Together they cover ions from H to U, all eight stopping units and energies from 10 eV (`LiSteelList`) to 10 GeV. For the other two inputs `SRModule.exe` stops with a run-time error and writes no table. In `AuH2O` (Au in water from 1 eV, correction 0.92) the straggling variance becomes negative, and in `HAu8G` (H in Au at 8 GeV) the exponential of Eq. (9) of [physics.md](physics.md) overflows. lunasrim raises `lunasrim.sr.SrError` for both.

## Running the tests

```sh
python -m unittest discover tests        # or: pytest
```

The suite runs lunasrim on all reference cases and checks the web application, the command line and the analytic energy straggling. On an Intel Xeon processor at 2.80 GHz the whole suite (125 tests) took 190 s.

## Processor dependence of the x87 arithmetic

The SRIM programs evaluate `Exp`, `Log`, `^`, `Atn`, `Sin`, `Cos` and `Sqr` with x87 instructions of the Visual Basic runtime, which runs the FPU with control word 0x137F (64-bit precision). We read this control word with winedbg inside TRIM's collision loop and in `SRModule.exe`. The results differ from the correctly rounded functions of Python's `math` module in the last bit for about 0.1 % of the arguments. In gas layers the deflections are small, so that $(1 - c_x^2)^{1/2}$ cancels and amplifies such a bit until a trajectory changes. Before we reproduced this arithmetic, 8 of the TRIM cases had last-digit differences, and in `HeArGasW` one ion took a visibly different path.

lunasrim runs the same x87 instruction sequences natively on x86-64 (`lunasrim/vbmath.py`, assembly source in `tools/vbmath.S`). We checked them against `msvbvm50.dll` itself: a small Windows program called the runtime's exports under wine with 600,000 random arguments, and lunasrim's routines agreed on all of them.

The x87 transcendental instructions (`fyl2x`, `f2xm1`, `fpatan`, `fsin`, `fcos`) are not correctly rounded, and their last bit differs between processor vendors. The reference outputs were made on an Intel processor. When `lunasrim.vbmath` is imported, it compares the native results for four arguments with the values that `msvbvm50.dll` gives on Intel processors, and it uses the native routines only if all four agree (`lunasrim.vbmath.NATIVE` is then true). On other processors, for example from AMD, on other architectures, or when executable memory cannot be allocated, it uses Python's `math` module. The square root is exactly specified, so it is run natively whenever possible and otherwise computed by an exact model of the 64-bit rounding. The environment variable `LUNASRIM_NO_NATIVE_MATH=1` forces the pure Python version.

With the pure Python math on Linux on this Intel machine, 45 of the 46 TRIM cases remain identical. In `HPbLong` one value of `LATERAL.txt` differs in its last digits; `tests/test_trim.py` allows this one line when `vbmath.NATIVE` is false. The result of the fallback also depends on the C library behind Python's `math` module. The GitHub runners of `.github/workflows/tests.yml` do not reproduce the Intel x87 values, so lunasrim uses the fallback there. On Linux (Python 3.9, 3.10, 3.12 and 3.13) and macOS (arm64, Python 3.12) all 46 TRIM cases were identical, including `HPbLong`. On Windows (Python 3.12) one more case differs: in `XeSiHigh/LATERAL.txt` two lateral spread values differ in the last digit (`12823989.E-04` against `12823988.E-04`), because the Microsoft C runtime rounds some transcendental results differently from glibc and the macOS library. `tests/test_trim.py` allows these two lines on Windows without native math. The SR module tables are identical on all platforms and with both kinds of math. The web application in the browser (Pyodide) has no native code either; there, 45 of the 46 TRIM cases are identical and `HPbLong` differs in one line of `LATERAL.txt`.

## Edge-case campaign

In addition to the kept references, we ran 57 TRIM inputs and 24 SR inputs once against the binaries. They are not kept in the repository. They covered:

- ion and energy extremes: 10 eV H, 1 GeV U, 89° incidence;
- targets with 8 layers, 1 Å layers, gas, solid and gas stacks, and compounds of 10 elements;
- every calculation type and every output-file switch;
- spellings of numbers in TRIM.IN (`2000`, `2000.0`, `2E3`, widths above 32767), which change Visual Basic's variable types;
- SR units, energy limits, compound corrections and list modes.

The campaign found two differences in the output files, which we fixed: the content of `TRIMOUT.txt` when only some per-ion files are switched on, and the header of the SR table for a compound correction of 0. After the fixes all cases were identical.

## Known limits

- Above $E/M_1 = 5.92$ GeV/u for protons and 5.85 GeV/u for heavier ions, the exponential of SRIM's phase correction overflows. `SRModule.exe` and `TRIM.exe` stop with Visual Basic run-time error 6. `lunasrim.sr.run` raises `SrError`, a TRIM calculation raises `lunasrim.stopping.VBOverflow`, and `lunasrim.energy_loss` raises `lunasrim.straggling.StragglingError`. `lunasrim.stopping_table` and the web application end the table at the last grid energy below the limit instead.
- On SRIM's energy grid the accumulated energy can miss E-Max by about $10^{-12}$ keV; `SRModule.exe` then integrates backwards and stops with a run-time error. `lunasrim.sr.run` raises `SrError` in that case, while `lunasrim.stopping_table` ends the table at E-Max.
- TRIM's special calculation types 4 to 8 (ion-only and neutron calculations) are partly present in the code but untested. Auto-save and the restart of an interrupted calculation are not implemented.
- The comparison covers the output files. TRIM's screen plots and the interactive SRIM programs other than TRIM and the SR module are not reproduced.

## Energy straggling

The analytic energy loss and straggling of [physics.md](physics.md), Section 6, is not part of SRIM, so it cannot be compared with SRIM output. We compare it instead with lunasrim's TRIM. For five ion-layer pairs the exit energies agree within 1.1 keV and the TRIM widths of 4000 transmitted ions are 1.5 to 7.0 % larger than the analytic widths; `tests/test_straggling.py` checks one of these cases with 600 ions. The table is in [physics.md](physics.md), Section 6.
