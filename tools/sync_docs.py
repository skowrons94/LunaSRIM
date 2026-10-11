#!/usr/bin/env python3
"""Copy docs/physics.md into the web application.

usage: tools/sync_docs.py [--check]

docs/physics.md is the single source of the physics documentation (GitHub
renders it, the web application shows it on its Physics tab).  The web
application reads lunasrim/web/static/physics.md, a committed copy, so that
the pip-installed package works without the repository.  Run this script
after editing docs/physics.md; tests/test_docs_page.py fails while the two
files differ.  With --check nothing is written and the exit status is 1 if
the copy is out of date.

Images referenced by the document (docs/img/...) are not copied into the
package: the page loads them from img/ next to index.html (the static site
built by tools/build_static.py has them) and otherwise from GitHub.
"""
import os
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, 'docs', 'physics.md')
COPY = os.path.join(ROOT, 'lunasrim', 'web', 'static', 'physics.md')


def images(path=DOC):
    """files of docs/img/ referenced by the Markdown document"""
    text = open(path, encoding='utf-8').read()
    refs = re.findall(r'!\[[^\]]*\]\(\s*<?([^)\s>]+)', text)
    refs += re.findall(r'<img[^>]+src=["\']([^"\']+)', text)
    out = []
    for r in refs:
        if '://' in r or r.startswith('/'):
            continue
        p = os.path.normpath(os.path.join(os.path.dirname(path), r.split('#')[0].split('?')[0]))
        if os.path.dirname(p) == os.path.join(ROOT, 'docs', 'img') and os.path.isfile(p) and p not in out:
            out.append(p)
    return out


def up_to_date():
    return (os.path.exists(COPY)
            and open(DOC, 'rb').read() == open(COPY, 'rb').read())


def sync():
    if not up_to_date():
        shutil.copyfile(DOC, COPY)
        print('copied %s -> %s' % (os.path.relpath(DOC, ROOT), os.path.relpath(COPY, ROOT)))
    else:
        print('%s is up to date' % os.path.relpath(COPY, ROOT))


def main(argv):
    if not os.path.exists(DOC):
        print('%s does not exist' % DOC, file=sys.stderr)
        return 1
    if '--check' in argv:
        ok = up_to_date()
        print('up to date' if ok else '%s differs from docs/physics.md: run tools/sync_docs.py'
              % os.path.relpath(COPY, ROOT))
        return 0 if ok else 1
    sync()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
