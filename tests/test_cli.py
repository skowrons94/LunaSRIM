"""Command line interface (lunasrim / python -m lunasrim)"""
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, ROOT)

import lunasrim                               # noqa: E402
from lunasrim import sr, straggling           # noqa: E402
from lunasrim.__main__ import main            # noqa: E402

REF_SR = os.path.join(HERE, 'reference_sr')
REF_TRIM = os.path.join(HERE, 'reference')


def run(*argv):
    """(exit code, stdout, stderr) of main(argv)"""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = main(list(argv))
        except SystemExit as ex:          # argparse
            code = ex.code
    return code, out.getvalue(), err.getvalue()


def _lines(path, skip):
    with open(path, 'rb') as fh:
        return [l for l in fh.read().decode('latin-1').split('\r\n') if skip not in l]


class TestCli(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='lunasrim_cli_')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_version_and_help(self):
        code, out, err = run('--version')
        self.assertEqual((code, out.strip()), (0, 'lunasrim %s' % lunasrim.__version__))
        for cmd in ('run', 'stopping', 'eloss', 'trim', 'web'):
            code, out, err = run(cmd, '--help')
            self.assertEqual(code, 0)
            self.assertIn('usage: lunasrim %s' % cmd, out)
        code, out, err = run()
        self.assertEqual(code, 2)
        self.assertIn('stopping', err)

    def test_module_entry_point(self):
        p = subprocess.run([sys.executable, '-m', 'lunasrim', '--version'], cwd=ROOT,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.strip(), 'lunasrim %s' % lunasrim.__version__)

    def test_sr_in_backwards_compatible(self):
        # "lunasrim SR.IN" writes the table named in SR.IN next to it, as SRModule.exe
        for case in ('HC', 'ArAirGas'):
            d = os.path.join(self.dir, case)
            shutil.copytree(os.path.join(REF_SR, case), d)
            si = sr.SrInput.read(os.path.join(d, 'SR.IN'))
            os.rename(os.path.join(d, si.output_name), os.path.join(d, 'reference.txt'))
            code, out, err = run(os.path.join(d, 'SR.IN'))
            self.assertEqual(code, 0, err)
            self.assertIn('stopping table written to', out)
            self.assertEqual(_lines(os.path.join(d, si.output_name), 'Calc. date'),
                             _lines(os.path.join(d, 'reference.txt'), 'Calc. date'))
        # -o puts it elsewhere; "run" is the same command
        code, out, err = run('run', os.path.join(REF_SR, 'HC', 'SR.IN'), '-o', self.dir, '-q')
        self.assertEqual((code, out), (0, ''))
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'HC.txt')))

    def test_trim_in_backwards_compatible(self):
        ref = os.path.join(REF_TRIM, 'ColFC1')
        outdir = os.path.join(self.dir, 'out')
        code, out, err = run(os.path.join(ref, 'TRIM.IN'), '-o', outdir)
        self.assertEqual(code, 0, err)
        self.assertIn('ions: 5', out)
        for name in ('RANGE.txt', 'VACANCY.txt', 'COLLISON.txt'):
            self.assertEqual(_lines(os.path.join(outdir, name), 'Calculation Date='),
                             _lines(os.path.join(ref, name), 'Calculation Date='), name)
        # -n before the file, old style
        code, out, err = run('-n', '2', '-q', os.path.join(ref, 'TRIM.IN'), '-o', outdir)
        self.assertEqual((code, out), (0, ''))

    def test_run_errors(self):
        code, out, err = run('no_such_file.in')
        self.assertEqual(code, 2)
        self.assertIn('no such file', err)
        code, out, err = run('run', 'a.in', 'b.in', '-o', self.dir)
        self.assertEqual(code, 2)

    def test_stopping(self):
        code, out, err = run('stopping', 'H', 'Si')
        self.assertEqual(code, 0, err)
        self.assertIn('MeV/(mg/cm2)', out)
        self.assertIn('10 keV', out)
        self.assertIn('1342 A', out)                          # range of 10 keV H in Si
        code, out, err = run('stopping', 'He', 'SiO2', '--density', '2.32', '--emin', '100', '--emax', '5MeV',
                             '--units', 'eV/A', '--csv')
        self.assertEqual(code, 0, err)
        lines = out.strip().split('\n')
        self.assertTrue(lines[0].startswith('energy_keV,S_elec [eV/A]'))
        rows = lunasrim.stopping_table('He', 'SiO2', density=2.32, emin=100, emax=5000, units=1)
        self.assertEqual(len(lines), len(rows) + 1)
        first = [float(x) for x in lines[1].split(',')]
        self.assertAlmostEqual(first[1] / rows[0][1], 1, places=6)
        # SRModule.exe format, to a file
        path = os.path.join(self.dir, 'HSi.txt')
        code, out, err = run('stopping', 'H', 'Si', '--srim', '-o', path)
        self.assertEqual(code, 0, err)
        with open(path, 'rb') as fh:
            text = fh.read().decode('latin-1')
        self.assertIn('Ion = Hydrogen', text)
        self.assertIn('\r\n', text)
        code, out, err = run('stopping', 'H', 'Si', '--srim', '--emin', '1MeV', '--emax', '1MeV')
        self.assertIn('Stopping Units =  MeV/(mg/cm2)', out)

    def test_stopping_errors(self):
        for args in (('stopping', 'H', 'SiO2'), ('stopping', 'Xx', 'Si'), ('stopping', 'H', 'Si', '--units', 'x'),
                     ('stopping', 'H', 'Si', '--emin', '1MeV', '--emax', '10keV')):
            code, out, err = run(*args)
            self.assertEqual(code, 2, args)
            self.assertTrue(err.startswith('lunasrim: error:'), err)
            self.assertNotIn('Traceback', err)

    def test_eloss(self):
        code, out, err = run('eloss', 'H', 'Si', '--energy', '400keV', '--thickness', '1um')
        self.assertEqual(code, 0, err)
        self.assertIn('330.464 keV', out)
        self.assertIn('4.4826 keV', out)
        code, out, err = run('eloss', 'H', 'Si', '--energy', '0.4MeV', '--thickness', '10000A', '--json')
        r = json.loads(out)
        r0 = straggling.energy_loss(1, 1.008, 400, [(14, 1, 28.086)], 2.3212, 1e4, 'A')
        self.assertEqual(r['energy_out'], r0['energy_out'])
        self.assertEqual(r['sigma'], r0['sigma'])
        self.assertNotIn('profile', r)
        code, out, err = run('eloss', 'He', 'N2', '--energy', '300', '--thickness', '20e15/cm2', '--json')
        r = json.loads(out)
        self.assertTrue(r['target']['gas'])
        self.assertAlmostEqual(r['areal_density'] / 2e16, 1, places=9)
        code, out, err = run('eloss', 'H', 'Si', '--energy', '400', '--thickness', '1')
        self.assertEqual(code, 2)
        self.assertIn('needs a unit', err)
        code, out, err = run('eloss', 'H', 'Si', '--energy', '10', '--thickness', '1mm')
        self.assertEqual(code, 0)
        self.assertIn('stop in the layer', out)

    def test_trim(self):
        outdir = os.path.join(self.dir, 'trim')
        code, out, err = run('trim', 'Ar', 'Si', '--energy', '50keV', '--thickness', '200nm', '--ions', '20',
                             '--seed', '7', '--out', outdir, '-q')
        self.assertEqual(code, 0, err)
        self.assertIn('mean projected range', out)
        self.assertIn('vacancies per ion', out)
        self.assertTrue(os.path.isfile(os.path.join(outdir, 'TRIM.IN')))
        self.assertTrue(os.path.isfile(os.path.join(outdir, 'SRIM Outputs', 'RANGE.txt')))
        res = lunasrim.trim('Ar', 'Si', energy=50, thickness=200, unit='nm', ions=20, seed=7,
                            directory=os.path.join(self.dir, 'api'))
        self.assertEqual(lunasrim.Results(os.path.join(outdir, 'SRIM Outputs')).range.average_range,
                         res.range.average_range)
        code, out, err = run('trim', 'Ar', 'Si', '--energy', '50keV', '--thickness', '200nm', '--calculation', '4')
        self.assertEqual(code, 2)


if __name__ == '__main__':
    unittest.main()
