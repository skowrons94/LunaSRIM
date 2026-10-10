"""Visual Basic runtime math (lunasrim.vbmath) against values computed by
msvbvm50.dll (SRIM-2013 setup) under wine, and native routines against the
pure Python fallback."""
import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from lunasrim import vbmath     # noqa: E402

# (function, arguments, value returned by msvbvm50.dll)
MSVBVM50 = [
    ('pow', (6432.957985524626, 2.0), 41382948.44352505),       # math.pow: ...046
    ('pow', (145265.48020290735, 2.0), 21102059738.581264),
    ('pow', (2.6329444762163088e-05, 2.0), 6.932396614837973e-10),
    ('atn', (0.09014091311591388,), 0.08989795332127876),         # math.atan: ...874
    ('pow', (2.0, 10.0), 1024.0),
    ('exp', (1.0,), 2.718281828459045),
    ('log', (10.0,), 2.302585092994046),
]


class TestVBMath(unittest.TestCase):

    @unittest.skipUnless(vbmath.NATIVE, 'native x87 routines not available')
    def test_msvbvm50_values(self):
        for name, args, v in MSVBVM50:
            self.assertEqual(getattr(vbmath, name)(*args), v, (name, args))

    def test_sqrt_64bit_precision(self):
        # TRIM runs with x87 control word 0x137f: fsqrt rounds to 64 bits first
        self.assertEqual(vbmath.sqrt(1 - 2 ** -53), 1.0)
        self.assertEqual(math.sqrt(1 - 2 ** -53), 1 - 2 ** -53)
        self.assertEqual(vbmath._sqrt64(1 - 2 ** -53), 1.0)

    def test_sqrt_native_vs_python(self):
        rnd = random.Random(3)
        xs = [1 - k * 2 ** -53 for k in range(1, 300)] + [1 + k * 2 ** -52 for k in range(1, 300)]
        xs += [rnd.uniform(0.25, 4) for i in range(3000)] + [10 ** rnd.uniform(-100, 100) for i in range(3000)]
        for x in xs:
            self.assertEqual(vbmath.sqrt(x), vbmath._sqrt64(x), x)

    def test_pow_special(self):
        self.assertEqual(vbmath.pow(-2.0, 3.0), -8.0)
        self.assertEqual(vbmath.pow(-2.0, 2.0), 4.0)
        self.assertEqual(vbmath.pow(0.0, 2.0), 0.0)
        self.assertEqual(vbmath.pow(5.0, 0.0), 1.0)
        self.assertEqual(vbmath.pow(10.0, -3), 0.001)


if __name__ == '__main__':
    unittest.main()
