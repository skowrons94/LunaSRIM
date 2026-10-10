"""Command line interface: ``python -m lunasrim TRIM.IN [-o DIR] [-n N]``"""
import argparse
import os
import sys
import time

from .trimin import TrimInput
from .srim import run_trim


def main(argv=None):
    p = argparse.ArgumentParser(
        prog='lunasrim',
        description='Run a SRIM-2013 TRIM calculation (TRIM.IN file) with the Python engine.')
    p.add_argument('trim_in', help='TRIM.IN input file')
    p.add_argument('-o', '--output', default='SRIM Outputs',
                   help='directory for the output files (default: "SRIM Outputs")')
    p.add_argument('-n', '--number-ions', type=int, default=None,
                   help='number of ions (default: the number given in TRIM.IN)')
    p.add_argument('-q', '--quiet', action='store_true', help='no progress output')
    a = p.parse_args(argv)
    ti = TrimInput.read(a.trim_in)
    n = int(ti.number_ions if a.number_ions is None else a.number_ions)
    t0 = time.time()
    state = {'last': 0.0}

    def progress(i, eng):
        now = time.time()
        if now - state['last'] > 1 or i == n:
            state['last'] = now
            sys.stderr.write('\rion %d / %d  (%.1f s)' % (i, n, now - t0))
            sys.stderr.flush()

    eng = run_trim(ti, a.output, number_ions=n, progress=None if a.quiet else progress)
    if not a.quiet:
        tal = eng.tal
        sys.stderr.write('\n')
        print('ions: %d  backscattered: %d  transmitted: %d' % (eng.ion_number, tal.n_back, tal.n_trans))
        print('average range %.6g A, straggling %.6g A' % (tal.g_a38, tal.g_a48))
        print('output written to %s' % os.path.abspath(a.output))
    return 0


if __name__ == '__main__':
    sys.exit(main())
