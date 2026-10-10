"""Compare lunasrim with output files of the original SRIM-2013 TRIM.exe.

Run with ``python -m unittest discover tests`` or ``pytest``.
"""
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))

from lunasrim.trimin import TrimInput     # noqa: E402
from lunasrim.srim import run_trim        # noqa: E402

REF = os.path.join(HERE, 'reference')

# lines known to differ from TRIM.exe: {case: {file: number of lines}}
KNOWN = {'HeAuC': {'TRANSMIT.txt': 2, 'TRIMOUT.txt': 2}}


def _lines(path):
    with open(path, 'rb') as fh:
        text = fh.read().decode('latin-1')
    return [l for l in text.split('\r\n') if 'Calculation Date=' not in l]


class TestAgainstTrim(unittest.TestCase):
    maxDiff = None

    def check_case(self, case):
        ref = os.path.join(REF, case)
        out = tempfile.mkdtemp(prefix='lunasrim_test_')
        try:
            run_trim(TrimInput.read(os.path.join(ref, 'TRIM.IN')), out)
            names = sorted(n for n in os.listdir(ref) if n.endswith('.txt'))
            self.assertTrue(names)
            for name in names:
                a = _lines(os.path.join(out, name))
                b = _lines(os.path.join(ref, name))
                self.assertEqual(len(a), len(b), '%s/%s: number of lines' % (case, name))
                bad = [(x, y) for x, y in zip(a, b) if x != y]
                allowed = KNOWN.get(case, {}).get(name, 0)
                self.assertLessEqual(len(bad), allowed, '%s/%s differs:\n%s' % (
                    case, name, '\n'.join('  lunasrim: %s\n  TRIM    : %s' % d for d in bad[:5])))
        finally:
            shutil.rmtree(out, ignore_errors=True)


def _make(case):
    return lambda self: self.check_case(case)


for _case in sorted(os.listdir(REF)):
    if os.path.isdir(os.path.join(REF, _case)):
        setattr(TestAgainstTrim, 'test_' + _case, _make(_case))


if __name__ == '__main__':
    unittest.main()
