#!/usr/bin/env python3
"""Build the static LunaSRIM website (GitHub Pages): no server, the lunasrim
package runs in the visitor's browser with Pyodide.

usage: tools/build_static.py OUT_DIR [--pyodide DIR_OR_TARBALL]

OUT_DIR receives the web application (index.html, app.js, ...), the
package as lunasrim.zip, and the Pyodide runtime in OUT_DIR/pyodide/
(pyodide-core, downloaded from the Pyodide GitHub release unless a local
directory or tarball is given).  Serve OUT_DIR with any static web server,
e.g. ``python3 -m http.server -d OUT_DIR``.
"""
import io
import os
import shutil
import sys
import tarfile
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, 'lunasrim')
STATIC = os.path.join(PKG, 'web', 'static')
PYODIDE_VERSION = '0.26.4'
PYODIDE_URL = ('https://github.com/pyodide/pyodide/releases/download/%s/pyodide-core-%s.tar.bz2'
               % (PYODIDE_VERSION, PYODIDE_VERSION))
FLAG = '<script>window.LUNASRIM_STATIC = true;</script>\n'


def package_zip():
    """lunasrim (Python sources and data, without the web front end)"""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
        for dirpath, dirnames, filenames in os.walk(PKG):
            dirnames[:] = [d for d in dirnames if d not in ('__pycache__', 'static')]
            for f in sorted(filenames):
                if f.endswith(('.pyc', '.pyo')):
                    continue
                p = os.path.join(dirpath, f)
                z.write(p, os.path.relpath(p, ROOT))
    return buf.getvalue()


def install_pyodide(dest, source=None):
    os.makedirs(dest, exist_ok=True)
    if source and os.path.isdir(source):
        for f in os.listdir(source):
            shutil.copy(os.path.join(source, f), dest)
        return
    if source:
        data = open(source, 'rb').read()
    else:
        print('downloading', PYODIDE_URL)
        data = urllib.request.urlopen(PYODIDE_URL).read()
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:bz2') as t:
        for m in t.getmembers():
            if m.isfile():
                name = os.path.basename(m.name)
                with open(os.path.join(dest, name), 'wb') as fh:
                    fh.write(t.extractfile(m).read())


def build(out, pyodide=None):
    if os.path.exists(out):
        shutil.rmtree(out)
    shutil.copytree(STATIC, out)
    index = os.path.join(out, 'index.html')
    html = open(index, encoding='utf-8').read()
    html = html.replace('<script src="app.js"></script>', FLAG + '<script src="app.js"></script>', 1)
    with open(index, 'w', encoding='utf-8') as fh:
        fh.write(html)
    with open(os.path.join(out, 'lunasrim.zip'), 'wb') as fh:
        fh.write(package_zip())
    install_pyodide(os.path.join(out, 'pyodide'), pyodide)
    open(os.path.join(out, '.nojekyll'), 'w').close()
    size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(out) for f in fs)
    print('static site in %s (%.1f MB)' % (out, size / 1e6))


def main(argv):
    if not argv or argv[0].startswith('-'):
        print(__doc__)
        return 2
    pyodide = None
    if '--pyodide' in argv:
        pyodide = argv[argv.index('--pyodide') + 1]
    build(argv[0], pyodide)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
