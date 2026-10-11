"""SR module (stopping and range tables) against SRModule.exe output files."""
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))

from lunasrim import sr     # noqa: E402

REF = os.path.join(HERE, 'reference_sr')


def _lines(path):
    with open(path, 'rb') as fh:
        return [l for l in fh.read().decode('latin-1').split('\\r\\n') if 'Calc. date' not in l]


class TestSR(unittest.TestCase):
    maxDiff = None

    def check_case(self, case):
        d = os.path.join(REF, case)
        si = sr.SrInput.read(os.path.join(d, 'SR.IN'))
        out = tempfile.mkdtemp(prefix='lunasrim_sr_')
        try:
            path = sr.run(os.path.join(d, 'SR.IN'), out)
            self.assertEqual(_lines(path), _lines(os.path.join(d, si.output_name)))
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_negative_variance(self):
        # SRModule.exe stops with a run-time error (empty output file)
        with self.assertRaises(sr.SrError):
            sr.run(os.path.join(REF, 'AuH2O', 'SR.IN'), tempfile.mkdtemp())

    def test_exp_overflow(self):
        # above 5.92 GeV/u SRModule.exe stops with run-time error 6 (Overflow)
        with self.assertRaises(sr.SrError):
            sr.run(os.path.join(REF, 'HAu8G', 'SR.IN'), tempfile.mkdtemp())

    def test_stopping_table(self):
        rows = sr.stopping_table(1, 1.008, [(6, 1, 12.011)], 2.253, 10, 1000, units=7)
        e, se, sn, rp, lon, lat = rows[0]
        self.assertAlmostEqual(e, 10.0, places=9)   # energies accumulate as in SRIM
        self.assertTrue(rp > 0 and se > sn > 0)


def _make(case):
    return lambda self: self.check_case(case)


for _case in sorted(os.listdir(REF)):
    if _case not in ('AuH2O', 'HAu8G'):
        setattr(TestSR, 'test_' + _case, _make(_case))

if __name__ == '__main__':
    unittest.main()
