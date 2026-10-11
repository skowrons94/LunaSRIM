"""API of the static website (lunasrim.web.pyodide_api, run here with CPython)."""
import json
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from lunasrim.web import pyodide_api as api     # noqa: E402

SI = {'name': 'Si', 'elements': [{'z': 14}], 'density': 2.321, 'width': 2000}


class TestStaticApi(unittest.TestCase):
    def call(self, method, path, body=None):
        st, text = api.handle(method, path, json.dumps(body) if body is not None else '')
        return st, json.loads(text)

    def test_data_and_stopping(self):
        st, els = self.call('GET', '/api/elements')
        self.assertEqual(st, 200)
        self.assertEqual(len(els), 92)
        st, info = self.call('GET', '/api/info')
        self.assertTrue(info['static'])
        st, d = self.call('POST', '/api/stopping', {'ion': {'z': 1}, 'target': SI, 'emin': 10, 'emax': 1000})
        self.assertEqual(st, 200)
        self.assertIn('Stopping Units', d['sr_output'])

    def test_errors(self):
        st, d = self.call('POST', '/api/stopping', {'ion': {'z': 1}, 'target': SI, 'emin': 0, 'emax': 1})
        self.assertEqual(st, 400)
        st, d = self.call('POST', '/api/stopping', {'ion': {'z': 99}, 'target': SI})
        self.assertEqual(st, 422)
        st, d = self.call('GET', '/api/nothing')
        self.assertEqual(st, 404)

    def test_empty_fields(self):
        # an empty or invalid input field arrives as null: a message, not a crash
        ion = {'z': 18, 'energy': None}
        layer = dict(SI, width=2000)
        for path, body, word in [
                ('/api/trim', {'ion': ion, 'layers': [layer]}, 'ion energy'),
                ('/api/trim', {'ion': dict(ion, energy=50), 'layers': [dict(layer, width=None)]}, 'layer width'),
                ('/api/stopping', {'ion': {'z': 1}, 'target': SI, 'emin': None, 'emax': 10}, 'minimum energy'),
                ('/api/implant', {'ion': {'z': 5}, 'target': SI, 'energies': [None]}, 'ion energy')]:
            st, d = self.call('POST', path, body)
            self.assertEqual(st, 422, (path, d))
            self.assertIn(word, d['detail'])
        # optional settings left empty take their defaults
        st, d = self.call('POST', '/api/trim', {'ion': dict(ion, energy=50), 'layers': [layer], 'number_ions': 5,
                                                'seed': None, 'plot_xmin': None, 'plot_xmax': None,
                                                'exyz': None, 'trajectories': None})
        self.assertEqual(st, 200, d)

    def test_trim_in_slices(self):
        st, j = self.call('POST', '/api/trim', {'ion': {'z': 18, 'energy': 50}, 'layers': [SI],
                                                 'number_ions': 20, 'seed': 5, 'trajectories': 3})
        self.assertEqual(st, 200)
        jid = j['id']
        n = 0
        while api.trim_step(jid, 0.01):
            n += 1
            st, s = self.call('GET', '/api/trim/%s?live=1' % jid)
            self.assertEqual(s['status'], 'running')
        st, s = self.call('GET', '/api/trim/' + jid)
        self.assertEqual(s['status'], 'finished')
        self.assertEqual(s['result']['ions'], 20)
        self.assertIn('RANGE.txt', s['result']['files'])
        self.assertTrue(api.trim_file(jid, 'TRIM.IN').startswith('LunaSRIM'))
        self.assertGreater(len(api.trim_zip(jid)), 1000)


if __name__ == '__main__':
    unittest.main()
