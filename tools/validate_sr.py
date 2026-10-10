#!/usr/bin/env python3
"""Compare lunasrim's SR module with SRModule.exe output files.

usage: tools/validate_sr.py CASE_DIR [CASE_DIR ...]

Each CASE_DIR holds an SR.IN and the file SRModule.exe wrote for it (see
tools/run_real_sr.sh).  The calculation date line is ignored.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from lunasrim import sr     # noqa: E402


def _lines(path):
    with open(path, 'rb') as fh:
        return [l for l in fh.read().decode('latin-1').split('\r\n') if 'Calc. date' not in l]


def compare(case):
    si = sr.SrInput.read(os.path.join(case, 'SR.IN'))
    ref = os.path.join(case, si.output_name)
    out = sr.run(os.path.join(case, 'SR.IN'), tempfile.mkdtemp(prefix='lunasrim_sr_'))
    a, b = _lines(out), _lines(ref)
    bad = [(i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return bad, len(a), len(b)


def main(argv):
    status = 0
    for case in argv:
        bad, na, nb = compare(case)
        ok = not bad and na == nb
        print('== %s: %s' % (case, 'identical' if ok else '%d lines differ (%d/%d lines)' % (len(bad), na, nb)))
        for i, x, y in bad[:4]:
            print('   %4d lunasrim: %s' % (i + 1, x))
            print('        SRModule: %s' % y)
        status = status or (0 if ok else 1)
    return status


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
