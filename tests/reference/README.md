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

Known differences: in HeAuC two direction cosines of transmitted ions
(TRANSMIT.txt/TRIMOUT.txt) differ in the 7th digit.
