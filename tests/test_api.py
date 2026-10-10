"""pysrim compatible API and output readers"""
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))

from lunasrim import Ion, Layer, Target, TRIM, Results, TrimInput   # noqa: E402

REF = os.path.join(HERE, 'reference')


class TestApi(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='lunasrim_api_')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_trim_in_and_results(self):
        ion = Ion('Ar', energy=50e3)
        self.assertEqual(ion.mass, 39.962)
        layer = Layer({'Si': {'stoich': 1, 'E_d': 15, 'lattice': 2, 'surface': 4.7}}, density=2.321, width=2000)
        trim = TRIM(Target([layer]), ion, number_ions=300, calculation=1, random_seed=777, ranges=True,
                    backscattered=True, transmit=True, sputtered=True)
        ti = trim.trim_input()
        ref = TrimInput.read(os.path.join(REF, 'ArSi', 'TRIM.IN'))
        for attr in ('ion_z', 'ion_mass', 'energy', 'number_ions', 'cascade', 'seed1'):
            self.assertEqual(getattr(ti, attr), getattr(ref, attr), attr)
        self.assertEqual(ti.elements[0].mass, 28.086)
        res = trim.run(self.dir)
        self.assertTrue(os.path.exists(os.path.join(self.dir, 'TRIM.IN')))
        real = Results(os.path.join(REF, 'ArSi'))
        for attr in ('ioniz', 'phonons', 'e2recoil', 'vacancy', 'range', 'lateral'):
            a, b = getattr(res, attr), getattr(real, attr)
            self.assertEqual(len(a.depth), 100)
            for x, y in zip(a.data, b.data):
                self.assertEqual(list(x), list(y), attr)
        self.assertEqual(res.range.average_range, real.range.average_range)
        self.assertEqual(len(res.range3d), 300)
        self.assertEqual(res.vacancy.vacancies_per_ion, 747)

    def test_readers_on_trim_files(self):
        r = Results(os.path.join(REF, 'HeAuC'))
        self.assertEqual(len(r.transmit), 298)
        self.assertTrue(all(z == 2 for z in r.transmit.z))
        r = Results(os.path.join(REF, 'NTaNFC'))
        self.assertEqual(len(r.vacancy.vacancies), 2)       # Ta and N columns
        self.assertGreater(len(r.sputter), 0)
        self.assertEqual(len(r.novac.number), 100)

    def test_exyz_and_collision_readers(self):
        r = Results(os.path.join(REF, 'ExyzBack'))
        e, x, y, z, se, t = r.exyz.ion(1)
        self.assertEqual(e[0], 2.0)                          # starts at E0 (keV)
        self.assertEqual(x[0], 0.0)
        r = Results(os.path.join(REF, 'ColGas2L'))
        self.assertTrue(any(v < 0 for v in r.exyz.positions[0]))   # backscattered ion
        r = Results(os.path.join(REF, 'ColFC1'))
        row = r.collisions.rows[0]
        self.assertEqual((row['ion'], row['atom'], row['displacements'], row['vacancies']), (1, 'Si', 2.0, 1.0))
        r = Results(os.path.join(REF, 'ColKP'))
        self.assertIn('vacancies', r.collisions.rows[0])

    def test_exyz_and_collisions_api(self):
        trim = TRIM(Target([Layer({'Si': {'stoich': 1}}, 2.321, 2000)]), Ion('Ar', 20e3), number_ions=3,
                    calculation=2, random_seed=5, collisions=2, exyz=100)
        res = trim.run(self.dir)
        self.assertGreater(len(res.exyz), 3)
        self.assertGreater(len(res.collisions), 0)

    def test_copy_output_files(self):
        trim = TRIM(Target([Layer({'C': {'stoich': 1}}, 2.0, 1000)]), Ion('H', 1e5), number_ions=5,
                    random_seed=1, transmit=True)
        trim.run(self.dir)
        dst = os.path.join(self.dir, 'copy')
        TRIM.copy_output_files(self.dir, dst)
        self.assertTrue(os.path.exists(os.path.join(dst, 'TRANSMIT.txt')))
        self.assertTrue(os.path.exists(os.path.join(dst, 'TRIM.IN')))


if __name__ == '__main__':
    unittest.main()
