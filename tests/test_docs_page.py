"""The Physics page of the web application shows docs/physics.md: the copy in
the package (lunasrim/web/static/physics.md, used by the pip-installed web
application) must be identical to it.  Run tools/sync_docs.py after editing
docs/physics.md."""
import os
import sys
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, ROOT)
DOC = os.path.join(ROOT, 'docs', 'physics.md')
COPY = os.path.join(ROOT, 'lunasrim', 'web', 'static', 'physics.md')
STATIC = os.path.join(ROOT, 'lunasrim', 'web', 'static')

try:
    from fastapi.testclient import TestClient
    from lunasrim.web.server import app
except ImportError:          # pragma: no cover
    app = None


@unittest.skipUnless(os.path.exists(DOC), 'docs/physics.md does not exist yet')
class TestPhysicsDoc(unittest.TestCase):
    def test_copy_is_identical(self):
        self.assertTrue(os.path.exists(COPY), 'lunasrim/web/static/physics.md missing: run tools/sync_docs.py')
        with open(DOC, 'rb') as a, open(COPY, 'rb') as b:
            self.assertTrue(a.read() == b.read(),
                            'lunasrim/web/static/physics.md differs from docs/physics.md: run tools/sync_docs.py')

    @unittest.skipIf(app is None, 'fastapi not installed')
    def test_server_route(self):
        r = TestClient(app).get('/physics.md')
        self.assertEqual(r.status_code, 200)
        self.assertTrue(r.headers['content-type'].startswith('text/markdown'))
        with open(DOC, encoding='utf-8') as fh:
            self.assertEqual(r.text, fh.read())


class TestPhysicsTab(unittest.TestCase):
    def test_tab_and_loader(self):
        with open(os.path.join(STATIC, 'index.html'), encoding='utf-8') as fh:
            html = fh.read()
        self.assertIn('data-view="physics"', html)
        self.assertIn('id="view-physics"', html)
        # the tab sits between Compare and About
        self.assertLess(html.index('data-view="compare"'), html.index('data-view="physics"'))
        self.assertLess(html.index('data-view="physics"'), html.index('data-view="about"'))
        with open(os.path.join(STATIC, 'app.js'), encoding='utf-8') as fh:
            js = fh.read()
        for s in ("fetch('physics.md'", 'vendor/marked.min.js', 'vendor/katex.min.js', 'vendor/katex.min.css',
                  'marked@12.0.2', 'katex@0.16.11'):
            self.assertIn(s, js)


if __name__ == '__main__':
    unittest.main()
