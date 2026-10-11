"""Simple API (lunasrim.api): name resolution, units, and equality with the
low-level functions"""
import math
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))

import lunasrim                                   # noqa: E402
from lunasrim import api, sr, straggling          # noqa: E402
from lunasrim.srim import Ion, Layer, Target, TRIM   # noqa: E402


class TestMaterial(unittest.TestCase):
    def test_element(self):
        m = lunasrim.material('Si')
        self.assertEqual((m.elements, m.density, m.gas, m.correction, m.source),
                         ([(14, 1.0, 28.086)], 2.3212, False, 1.0, 'element'))
        for t in ('si', 'silicon', 14, '14'):
            self.assertEqual(lunasrim.material(t).elements, m.elements)
        self.assertEqual(lunasrim.material('Au').density, 19.311)

    def test_gas_element(self):
        for t in ('N', 'N2', 'nitrogen'):
            m = lunasrim.material(t)
            self.assertTrue(m.gas)
            self.assertEqual(m.density, 0.001251)          # SRIM's gas density
            self.assertEqual(m.elements[0][0], 7)
        self.assertEqual(lunasrim.material('N2').elements, [(7, 2.0, 14.007)])
        self.assertEqual(lunasrim.material('Ar').density, 0.001784)
        m = lunasrim.material('N', gas=False)
        self.assertFalse(m.gas)
        self.assertEqual(m.density, 1.026)                 # solid density
        with self.assertRaises(ValueError):
            lunasrim.material('Si', gas=True)              # no SRIM gas density for Si
        self.assertTrue(lunasrim.material('Si', gas=True, density=1e-3).gas)

    def test_formula(self):
        m = lunasrim.material('SiO2', density=2.32)
        self.assertEqual((m.elements, m.density, m.gas, m.source),
                         ([(14, 1.0, 28.086), (8, 2.0, 15.999)], 2.32, False, 'formula'))
        with self.assertRaises(ValueError) as cm:
            lunasrim.material('SiO2')
        self.assertIn('density', str(cm.exception))
        self.assertIn('quartz', str(cm.exception))          # suggests the compound
        self.assertEqual([z for z, s, m in lunasrim.material('C8H8', density=1.05).elements], [6, 1])
        self.assertEqual(lunasrim.material('CH3COOH', density=1.05).elements,
                         [(6, 2.0, 12.011), (1, 4.0, 1.008), (8, 2.0, 15.999)])
        # 'CO' is carbon monoxide, not cobalt
        self.assertEqual([z for z, s, m in lunasrim.material('CO', density=1.0).elements], [6, 8])
        self.assertEqual(lunasrim.material('Co').elements[0][0], 27)

    def test_compound(self):
        m = lunasrim.material('Kapton')
        self.assertEqual((m.name, m.density, m.gas, m.source),
                         ('Kapton Polyimide Film (ICRU-179)', 1.42, False, 'compound'))
        self.assertEqual([z for z, s, ms in m.elements], [1, 6, 7, 8])
        self.assertAlmostEqual(sum(s for z, s, ms in m.elements), 1.0)
        self.assertEqual(lunasrim.material('mylar').density, 1.397)
        self.assertEqual(lunasrim.material('Water_Liquid').density, 1.0)
        self.assertEqual(lunasrim.material('water_liquid (icru-276)').density, 1.0)
        air = lunasrim.material('Air')
        self.assertTrue(air.gas)
        self.assertEqual(air.density, 0.00120484)
        self.assertEqual(lunasrim.material('Al2O3').density, 3.98)        # dictionary entry
        self.assertEqual(lunasrim.material('Kapton', density=1.4, gas=True).density, 1.4)
        with self.assertRaises(ValueError) as cm:
            lunasrim.material('Water')                      # liquid or vapour?
        self.assertIn('Water_vapor', str(cm.exception))

    def test_dict_and_errors(self):
        m = lunasrim.material({'Si': 1, 'O': 2}, density=2.2)
        self.assertEqual(m.elements, [(14, 1.0, 28.086), (8, 2.0, 15.999)])
        self.assertIs(lunasrim.material(m).__class__, api.Material)
        for bad in ('Xx', 'Unobtainium', '', 'H2O-'):
            with self.assertRaises(ValueError):
                lunasrim.material(bad, density=1.0)
        with self.assertRaises(ValueError):
            lunasrim.material('Si', density=0)

    def test_ion(self):
        self.assertEqual(lunasrim.resolve_ion('He'), (2, 4.003))
        self.assertEqual(lunasrim.resolve_ion(18), (18, 39.962))
        self.assertEqual(lunasrim.resolve_ion('ar', mass=40), (18, 40.0))
        with self.assertRaises(ValueError):
            lunasrim.resolve_ion('Q')
        with self.assertRaises(ValueError):
            lunasrim.resolve_ion(93)


class TestUnits(unittest.TestCase):
    def test_stopping_units(self):
        names = {1: ['eV/A', 'eV/Angstrom', 'ev / a', 'eV/Å'], 2: ['keV/um', 'keV/micron', 'keV/µm'],
                 3: ['MeV/mm'], 4: ['keV/(ug/cm2)', 'keV/ug/cm2'], 5: ['MeV/(mg/cm2)'],
                 6: ['keV/(mg/cm2)'], 7: ['eV/(1E15 atoms/cm2)', 'eV/(1e15/cm2)'], 8: ['LSS', 'L.S.S.']}
        for u, ns in names.items():
            self.assertEqual(api.stopping_units(u), u)
            self.assertEqual(api.stopping_units(str(u)), u)
            for n in ns:
                self.assertEqual(api.stopping_units(n), u, n)
        for bad in (0, 9, 'eV/nm', 'x'):
            with self.assertRaises(ValueError):
                api.stopping_units(bad)

    def test_energy(self):
        for text, kev in (('400', 400), ('400keV', 400), ('2MeV', 2000), ('2 mev', 2000), ('50eV', 0.05),
                          ('1GeV', 1e6), (3, 3), ('1.5e3 keV', 1500)):
            self.assertEqual(api.parse_energy(text), kev)
        for bad in ('x', '4x', '4 keVs', ''):
            with self.assertRaises(ValueError):
                api.parse_energy(bad)

    def test_thickness(self):
        cases = {'1um': (1, 'um'), '200nm': (200, 'nm'), '50A': (50, 'A'), '100ug/cm2': (100, 'ug/cm2'),
                 '0.5mg/cm2': (0.5, 'mg/cm2'), '2e19atoms/cm2': (2e19, 'atoms/cm2'),
                 '20e15/cm2': (2e16, 'atoms/cm2'), '10 1e15/cm2': (10, '1e15/cm2'), '22cm': (220, 'mm'),
                 '3 mm': (3, 'mm'), '1 µm': (1, 'um')}
        for text, res in cases.items():
            self.assertEqual(api.parse_thickness(text), res, text)
        self.assertEqual(api.parse_thickness(5, 'nm'), (5, 'nm'))
        for bad in ('1', '1 parsec', 'um'):
            with self.assertRaises(ValueError):
                api.parse_thickness(bad)
        for u in straggling.THICKNESS_UNITS:
            self.assertEqual(api.parse_thickness('1 ' + u)[1], u)


class TestSameAsLowLevel(unittest.TestCase):
    def test_stopping_table(self):
        self.assertEqual(lunasrim.stopping_table('H', 'Si', emin=10, emax=10000),
                         sr.stopping_table(1, 1.008, [(14, 1, 28.086)], 2.3212, 10, 10000, units=5))
        self.assertEqual(lunasrim.stopping_table('He', 'SiO2', density=2.32, emin=100, emax=5000, units='eV/A'),
                         sr.stopping_table(2, 4.003, [(14, 1, 28.086), (8, 2, 15.999)], 2.32, 100, 5000, units=1))
        self.assertEqual(lunasrim.stopping_table('He', 'N2', emin='1MeV', emax='2MeV', units=7),
                         sr.stopping_table(2, 4.003, [(7, 2, 14.007)], 0.001251, 1000, 2000, phase=1, units=7))
        with self.assertRaises(ValueError):
            lunasrim.stopping_table('H', 'Si', emin=100, emax=10)

    def test_stopping_report_is_sr_module(self):
        # same table as the SRModule.exe reference case H in C (tests/reference_sr/HC)
        d = os.path.join(HERE, 'reference_sr', 'HC')
        si = sr.SrInput.read(os.path.join(d, 'SR.IN'))
        self.assertEqual([(e.z, e.stoich) for e in si.elements], [(6, 1)])
        t = api.stopping_report(si.ion_z, si.elements[0].z, emin=si.emin, emax=si.emax, units=si.units,
                                mass=si.ion_mass, density=si.density, gas=bool(si.phase))
        ref = sr.StoppingTable.from_file(os.path.join(d, 'SR.IN'))
        self.assertEqual(t.rows, ref.rows)

    def test_energy_loss(self):
        r = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um')
        r0 = straggling.energy_loss(1, 1.008, 400, [(14, 1, 28.086)], density=2.3212, thickness=1, unit='um')
        for k in r0:
            self.assertEqual(r[k], r0[k], k)
        self.assertAlmostEqual(r['energy_out'], 330.464, places=3)
        self.assertAlmostEqual(r['sigma'], 4.4826, places=4)
        self.assertEqual(lunasrim.energy_loss('H', 'Si', energy='0.4MeV', thickness='1um')['sigma'], r0['sigma'])
        g = lunasrim.energy_loss('He', 'N2', energy=300, thickness=2e17, unit='atoms/cm2', gas=True)
        g0 = straggling.energy_loss(2, 4.003, 300, [(7, 2, 14.007)], 0.001251, 2e17, 'atoms/cm2', phase=1)
        self.assertEqual(g['energy_out'], g0['energy_out'])
        self.assertEqual(g['sigma'], g0['sigma'])
        # tilted layer: path = thickness / cos(angle)
        t = lunasrim.energy_loss('H', 'Si', energy=400, thickness=1, unit='um', angle=60)
        t0 = straggling.energy_loss(1, 1.008, 400, [(14, 1, 28.086)], 2.3212, 1 / math.cos(math.radians(60)), 'um')
        self.assertEqual(t['energy_out'], t0['energy_out'])
        self.assertAlmostEqual(t['path_factor'], 2.0)


class TestTrim(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='lunasrim_api_simple_')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_trim_same_as_pysrim_api(self):
        res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=40, seed=4321,
                            directory=os.path.join(self.dir, 'a'))
        self.assertEqual(res.directory, os.path.join(self.dir, 'a', 'SRIM Outputs'))
        layer = Layer({'Si': {'stoich': 1.0, 'E_d': 15, 'lattice': 3.0, 'surface': 4.7}}, density=2.3212,
                      width=2000.0, name='Si')
        ref = TRIM(Target([layer]), Ion('Ar', energy=50e3), number_ions=40, calculation=1, random_seed=4321,
                   ranges=True, backscattered=True, transmit=True, sputtered=True,
                   description='Ar (50 keV) into Si').run(os.path.join(self.dir, 'b'))
        self.assertEqual(res.range.average_range, ref.range.average_range)
        self.assertEqual(res.range.straggling, ref.range.straggling)
        self.assertEqual(res.vacancy.vacancies_per_ion, ref.vacancy.vacancies_per_ion)
        self.assertGreater(res.range.average_range, 300)
        with open(os.path.join(self.dir, 'a', 'TRIM.IN')) as fh:
            self.assertIn('Si', fh.read())

    def test_trim_temporary_directory_and_errors(self):
        res = lunasrim.trim('H', 'Kapton', energy=10, thickness='100nm', ions=3, seed=1, calculation=2)
        try:
            self.assertTrue(os.path.isfile(os.path.join(res.directory, 'RANGE.txt')))
        finally:
            shutil.rmtree(os.path.dirname(res.directory), ignore_errors=True)
        with self.assertRaises(ValueError):
            lunasrim.trim('H', 'Si', energy=10, thickness=100, unit='nm', calculation=5)
        with self.assertRaises(ValueError):
            lunasrim.trim('H', 'Si', energy=10, thickness=0, unit='nm')


if __name__ == '__main__':
    unittest.main()
