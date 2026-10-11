"""Energy loss and straggling (lunasrim.straggling) against Bohr's formula
and against the TRIM Monte Carlo (transmitted ions)."""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from lunasrim.straggling import energy_loss, BOHR_K, StragglingError     # noqa: E402
from lunasrim.srim import Ion, Layer, Target, TRIM                         # noqa: E402

SI = [(14, 1, 28.086)]


class TestStraggling(unittest.TestCase):

    def test_thin_layer_is_bohr(self):
        # 10 nm: the stopping slope does not matter yet
        r = energy_loss(1, 1.008, 2000, SI, 2.3212, 10, 'nm')
        n = 0.6022 * 2.3212 / 28.086
        bohr = math.sqrt(BOHR_K * 14 * n * 100) / 1000
        self.assertAlmostEqual(r['sigma_bohr'] / bohr, 1, delta=0.01)
        self.assertAlmostEqual(r['sigma'] / bohr, 1, delta=0.01)
        self.assertAlmostEqual(r['fwhm'] / r['sigma'], 2.35482, places=4)

    def test_units(self):
        a = energy_loss(1, 1.008, 400, SI, 2.3212, 1, 'um')
        for v, u in [(1000, 'nm'), (232.12, 'ug/cm2'), (a['areal_density'] / 1e15, '1e15/cm2')]:
            b = energy_loss(1, 1.008, 400, SI, 2.3212, v, u)
            self.assertAlmostEqual(a['energy_out'], b['energy_out'], places=6)
            self.assertAlmostEqual(a['sigma'], b['sigma'], places=6)

    def test_stopping_power_integral(self):
        # thin layer: energy loss = S(E) * thickness
        r = energy_loss(2, 4.003, 2000, [(79, 1, 196.97)], 19.311, 5, 'nm')
        s = r['profile']['se'][0] + r['profile']['sn'][0]
        self.assertAlmostEqual(r['energy_loss'] / (s * 50 / 1000), 1, delta=0.01)

    def test_ions_stop(self):
        r = energy_loss(2, 4.003, 100, SI, 2.3212, 10, 'um')
        self.assertTrue(r['stopped'])
        self.assertLess(r['depth_reached'], 1e5)
        with self.assertRaises(StragglingError):
            energy_loss(1, 1.008, 100, SI, 2.3212, 1, 'furlong')

    def test_against_trim(self):
        # 400 keV protons through 1 um of silicon
        r = energy_loss(1, 1.008, 400, SI, 2.3212, 10000, 'A')
        lay = Layer({'Si': {'stoich': 1, 'E_d': 15.0, 'lattice': 2.0, 'surface': 4.7}}, density=2.3212, width=10000)
        res = TRIM(Target([lay]), Ion('H', energy=400e3, mass=1.008), number_ions=600, calculation=1,
                   transmit=True, random_seed=11).run(tempfile.mkdtemp())
        e = sorted(x / 1000 for x in res.transmit.energy)
        n = len(e)
        sigma = (e[int(0.8413 * n)] - e[int(0.1587 * n)]) / 2       # without the nuclear tail
        self.assertAlmostEqual(e[n // 2] / r['energy_out'], 1, delta=0.01)
        self.assertAlmostEqual(sigma / r['sigma'], 1, delta=0.15)


if __name__ == '__main__':
    unittest.main()
