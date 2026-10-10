#!/usr/bin/env python3
"""Compare lunasrim with the original TRIM.exe output files.

usage: tools/validate.py CASE_DIR [CASE_DIR ...]

Each CASE_DIR contains a TRIM.IN and the output files written by SRIM-2013
for it (see tools/run_real_trim.sh).  lunasrim is run on the same TRIM.IN
and every output file is compared line by line (the calculation date/time
line is ignored).  Exit status 1 if any file differs.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from lunasrim.trimin import TrimInput      # noqa: E402
from lunasrim.srim import run_trim         # noqa: E402

FILES = ('RANGE.txt', 'IONIZ.txt', 'PHONON.txt', 'E2RECOIL.txt', 'VACANCY.txt', 'NOVAC.txt', 'EXYZ.txt', 'COLLISON.txt',
         'LATERAL.txt', 'TDATA.txt', 'TRIMOUT.txt', 'TRANSMIT.txt', 'BACKSCAT.txt',
         'SPUTTER.txt', 'RANGE_3D.txt')


def _lines(path):
    with open(path, 'rb') as fh:
        return [l for l in fh.read().decode('latin-1').split('\r\n') if 'Calculation Date=' not in l and not l.startswith('Calc. Date=')]


def compare(case, verbose=True):
    ti = TrimInput.read(os.path.join(case, 'TRIM.IN'))
    out = tempfile.mkdtemp(prefix='lunasrim_')
    run_trim(ti, out)
    bad = 0
    report = []
    for name in FILES:
        ref = os.path.join(case, name)
        if not os.path.exists(ref):
            continue
        mine = os.path.join(out, name)
        if not os.path.exists(mine):
            report.append('%s: missing' % name)
            bad += 1
            continue
        a, b = _lines(mine), _lines(ref)
        nd = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
        report.append('%s: %s' % (name, 'ok' if nd == 0 else '%d lines differ' % nd))
        if nd:
            bad += 1
            if verbose:
                for x, y in zip(a, b):
                    if x != y:
                        report.append('    lunasrim: ' + x)
                        report.append('    TRIM    : ' + y)
                        break
    extra = sorted(n for n in FILES if os.path.exists(os.path.join(out, n)) and
                   not os.path.exists(os.path.join(case, n)))
    for name in extra:
        report.append('%s: written by lunasrim but not by TRIM' % name)
        bad += 1
    return bad, report, out


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    status = 0
    for case in argv:
        bad, report, out = compare(case)
        print('== %s (%s)' % (case, 'identical' if bad == 0 else '%d files differ' % bad))
        for l in report:
            print('  ' + l)
        status = status or (1 if bad else 0)
    return status


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
