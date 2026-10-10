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
| HeAuC  | Kinchin-Pease | two layers, transmission |
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

| AuCKP       | Kinchin-Pease | 2 MeV Au in graphite (heavy ion, light target) |
| HPbLong     | Kinchin-Pease | 10 MeV H through 0.4 mm Pb |
| NeSiGlance  | full cascade  | 1 keV Ne on Si at 85 deg (mostly backscattered) |
| HeBragg     | Kinchin-Pease | water / SiO2 with Bragg corrections 0.92 and 1.07 |
| SiSiO2Si    | full cascade  | Si / SiO2 / Si, the same element in three layers |
| ArWLow      | monolayer / sputtering | 100 eV Ar on W |
| XeSiHigh    | Kinchin-Pease | 200 MeV Xe in Si |
| LiNiWin     | Kinchin-Pease | plot window 5000 .. 17500.5 A, 30 deg |
| CSteel8     | Kinchin-Pease | 3 MeV C in an 8-element steel |
| UAgFCSP     | full cascade  | 50 keV U on Ag / AgAu, 30 deg |
| HeArGasW    | Kinchin-Pease | 2 MeV He through 1 mm Ar gas into W |
| NAlGaAsSP   | monolayer / sputtering | 5 keV N on AlGaAs at 60 deg |
| ColKP       | Kinchin-Pease | COLLISON.txt (ion collisions with Kinchin-Pease vacancies) |
| ColKPcomp   | Kinchin-Pease | COLLISON.txt in GaAs, 20 deg |
| ColFC1      | full cascade  | COLLISON.txt, collisions with damage counts |
| ColFC1comp  | full cascade  | COLLISON.txt in TaN, all per-ion files too |
| ColFC       | full cascade  | COLLISON.txt with the listing of every cascade |
| ColFC2Lay   | full cascade  | cascade listing and EXYZ.txt, SiO2 / Si |
| ColSP       | monolayer / sputtering | cascade listing, 1 keV Ar on Cu |
| ColSPcomp   | monolayer / sputtering | cascade listing, 2 keV Xe on InP |
| ColGas2L    | full cascade  | cascade listing, EXYZ.txt, TaN / N2 gas, 30 deg |
| ExyzKP      | Kinchin-Pease | EXYZ.txt only (no TRIMOUT.txt is written) |
| ExyzFCtr    | full cascade  | EXYZ.txt, transmitted ions |
| ExyzBack    | Kinchin-Pease | EXYZ.txt, backscattered ions (2 keV H on Au) |

All files agree line by line (lunasrim also writes no file that TRIM does
not write) with lunasrim's native VB runtime math
(`lunasrim/vbmath.py`).  Before it was used, HeAuC, UFeKP, HKaptonBR,
HeMylarAir, HeArGasW, HPbLong, CSteel8 and XeSiHigh had a few last-digit
differences, and in HeArGasW one ion took a visibly different path.  The
cause was found by running TRIM.exe under winedbg (ion 166 of HeArGasW,
reproduced as a one-ion run by giving its random number state as the two
seeds of TRIM.IN): the x87 control word is 0x137f (64-bit precision), so
`Sqr` rounds twice, and `Exp`, `^`, `Log` and `Atn` are x87 sequences of the
VB runtime.
