Output files of the original SRIM-2013 `TRIM.exe` (run under wine in the
`costrouc/pysrim` docker image with `tools/run_real_trim.sh`) for the
`TRIM.IN` in each directory.  `tests/test_trim.py` checks that lunasrim
reproduces them line by line (apart from the calculation date).

| case   | calculation | what it covers |
|--------|-------------|----------------|
| ArSi   | Kinchin-Pease | 300 Ar 50 keV in Si |
| AuSiKP | Kinchin-Pease | heavy ion, 1 MeV Au in Si |
| ArSiFC | full cascade  | recoil cascades |
| HW     | full cascade  | 45 deg incidence, backscattering |
| HeAuC  | Kinchin-Pease | two layers, transmission (only 2 last-digit differences, see below) |
| HeN    | Kinchin-Pease | gas target, 22 cm wide layer |
| NTaN   | full cascade  | compound + gas layer, 60 deg, sputtering |
| NTaNFC | full cascade  | as NTaN, full cascades |
| XeGaAs | Kinchin-Pease | four elements in two compounds |
| ArCuSP | monolayer / sputtering (calculation 3) | |
| BSiPW  | Kinchin-Pease | plot window (Xmin, Xmax) |
| OAlFC  | full cascade  | three layers with non integer widths |
| UFeKP  | Kinchin-Pease | 50 MeV U in steel (fractional stoichiometries) |
| HKaptonBR   | Kinchin-Pease | 500 keV H in Kapton (C22H10N2O5), transmission |
| OBoneFC     | full cascade  | O in bone, six elements, fractional stoichiometries |
| ArGaNSP     | monolayer / sputtering | 2 keV Ar on GaN at 45 deg |
| HeMylarAir  | Kinchin-Pease | 5 MeV He through Mylar / 0.5 mm air (gas) / Si, missing surface energies (defaults) |
| WebSiO2     | Kinchin-Pease | quartz from the compound dictionary (as built by the web app) |
| WebAirPMMA  | Kinchin-Pease | PMMA / dry air (gas) / GaAs, 7 elements, 10 deg |
| WebKaptonFC | full cascade  | Si in Kapton from the compound dictionary |

Known differences (last printed digit only):

* HeAuC: two direction cosines of transmitted ions (TRANSMIT.txt/TRIMOUT.txt).
* HKaptonBR, HeMylarAir: direction cosines of some transmitted ions
  (TRANSMIT.txt/TRIMOUT.txt, 7th decimal); positions and energies agree.
  In HeMylarAir the ions cross 0.5 mm of gas, so the round-off of the many
  direction updates accumulates.
* UFeKP: the lateral position of ion 29 (RANGE_3D.txt, 5th digit) and two
  LATERAL.txt values that depend on it.

For UFeKP the whole trajectory of ion 29 (3599 collisions) was checked
against TRIM.exe running in an x86 emulator with Python implementations of
the VB runtime math functions: they agree bit for bit.  The differences
with TRIM.exe under wine therefore come from the last bits of the
floating point library functions (x87 / msvbvm50 versus the C library used
by Python), which occasionally change a printed digit.
