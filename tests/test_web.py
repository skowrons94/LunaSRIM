"""Web application API (skipped when fastapi is not installed)."""
import os
import sys
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

try:
    from fastapi.testclient import TestClient
    from lunasrim.web.server import app
except ImportError:          # pragma: no cover
    app = None

from lunasrim import sr     # noqa: E402


@unittest.skipIf(app is None, 'fastapi not installed')
class TestWeb(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c = TestClient(app)

    def test_data(self):
        els = self.c.get('/api/elements').json()
        self.assertEqual(len(els), 92)
        self.assertEqual(els[13]['symbol'], 'Si')
        comps = self.c.get('/api/compounds').json()
        self.assertGreater(len(comps), 300)

    def test_stopping_matches_sr_module(self):
        q = {'ion': {'z': 1, 'mass': 1.008}, 'target': {'elements': [{'z': 6, 'mass': 12.011}], 'density': 2.253},
             'emin': 10, 'emax': 1000, 'units': 7}
        d = self.c.post('/api/stopping', json=q).json()
        ref = sr.stopping_table(1, 1.008, [(6, 1, 12.011)], 2.253, 10, 1000, units=7)
        self.assertEqual(len(d['rows']), len(ref))
        self.assertAlmostEqual(d['rows'][-1]['range'], ref[-1][3])
        self.assertIn('Calculation using SRIM-2006', d['sr_output'])

    def test_implant(self):
        q = {'ion': {'z': 5}, 'target': {'elements': [{'z': 14}], 'density': 2.321}, 'energies': [10, 30],
             'doses': [1e14, 1e15]}
        d = self.c.post('/api/implant', json=q).json()
        self.assertEqual(len(d['results']), 2)
        self.assertLess(d['results'][0]['range'], d['results'][1]['range'])
        self.assertEqual(len(d['depth']), len(d['total']))

    def test_bad_input(self):
        q = {'ion': {'z': 1}, 'target': {'elements': [{'z': 6}], 'density': 2.2}, 'emin': 10, 'emax': 1}
        self.assertEqual(self.c.post('/api/stopping', json=q).status_code, 400)

    def test_trim_job(self):
        q = {'ion': {'z': 18, 'energy': 50}, 'number_ions': 30, 'calculation': 1, 'seed': 777, 'trajectories': 5,
             'layers': [{'name': 'Si', 'elements': [{'z': 14, 'mass': 28.086, 'displacement': 15, 'lattice': 2,
                                                     'surface': 4.7}], 'density': 2.321, 'width': 2000}]}
        jid = self.c.post('/api/trim', json=q).json()['id']
        for i in range(600):
            s = self.c.get('/api/trim/' + jid).json()
            if s['status'] in ('finished', 'error', 'cancelled'):
                break
            time.sleep(0.1)
        self.assertEqual(s['status'], 'finished', s.get('error'))
        r = s['result']
        self.assertEqual(r['ions'], 30)
        self.assertEqual(len(r['trajectories']), 5)
        self.assertIn('RANGE.txt', r['files'])
        z = self.c.get('/api/trim/%s/files.zip' % jid)
        self.assertEqual(z.status_code, 200)

    def test_trim_exyz_collisions(self):
        q = {'ion': {'z': 2, 'energy': 500}, 'number_ions': 4, 'calculation': 2, 'seed': 3, 'trajectories': 0,
             'exyz': 1000, 'collisions': 2,
             'layers': [{'name': 'Si', 'elements': [{'z': 14}], 'density': 2.321, 'width': 20000}]}
        jid = self.c.post('/api/trim', json=q).json()['id']
        for i in range(600):
            s = self.c.get('/api/trim/' + jid).json()
            if s['status'] in ('finished', 'error', 'cancelled'):
                break
            time.sleep(0.1)
        self.assertEqual(s['status'], 'finished', s.get('error'))
        r = s['result']
        self.assertIn('EXYZ.txt', r['files'])
        self.assertIn('COLLISON.txt', r['files'])
        self.assertEqual(len(r['exyz']), 4)
        self.assertEqual(r['exyz'][0]['e'][0], 500.0)


if __name__ == '__main__':
    unittest.main()
