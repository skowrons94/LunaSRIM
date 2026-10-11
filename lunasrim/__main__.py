"""Command line interface (``lunasrim`` or ``python -m lunasrim``).

::

    lunasrim TRIM.IN [-o DIR] [-n N]        TRIM calculation of a TRIM.IN file
    lunasrim SR.IN [-o DIR]                 stopping / range table of an SR.IN file
    lunasrim run FILE [FILE ...]            the same, for one or more files
    lunasrim stopping ION TARGET            stopping and range table
    lunasrim eloss ION TARGET --energy E --thickness T   energy loss and straggling
    lunasrim trim ION TARGET --energy E --thickness T    TRIM Monte Carlo in one layer
    lunasrim web                            start the web application
    lunasrim --version

Energies are in keV, or with a unit (``400keV``, ``2MeV``); thicknesses
have a unit (``1um``, ``200nm``, ``50A``, ``100ug/cm2``, ``2e19atoms/cm2``).
"""
import argparse
import csv
import io
import json
import os
import sys
import time

from . import __version__
from . import api

COMMANDS = ('run', 'stopping', 'eloss', 'trim', 'web')


class UserError(Exception):
    pass


# ----------------------------------------------------------------- formatting
def _fmt_energy(kev):
    for unit, f in (('GeV', 1e6), ('MeV', 1e3), ('keV', 1.0)):
        if kev >= f * 0.999999:
            return '%.4g %s' % (kev / f, unit)
    return '%.4g eV' % (kev * 1e3)


def _fmt_length(a):
    for unit, f in (('m', 1e10), ('mm', 1e7), ('um', 1e4)):
        if abs(a) >= f:
            return '%.4g %s' % (a / f, unit)
    return '%.4g A' % a


def _ion_text(z, m):
    from . import data_files
    return '%s (Z=%d, %g amu)' % (data_files.atomdata(z).symbol, z, m)


def _material_text(mat):
    return '%s [%s], %.6g g/cm3, %s' % (mat.name, mat.formula(), mat.density, 'gas' if mat.gas else 'solid')


def _write(text, path):
    if path:
        with open(path, 'w', newline='') as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)


def _target_args(p):
    p.add_argument('ion', help='ion: element symbol, name or Z (e.g. H, He, 18)')
    p.add_argument('target', help='target: element (Si), formula (SiO2, needs --density) or SRIM compound '
                                  'name (Kapton, Mylar, Water_Liquid)')
    p.add_argument('--density', type=float, metavar='G_CM3',
                   help='target density in g/cm3 (default: SRIM value for elements and compounds)')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--gas', dest='gas', action='store_true', default=None,
                   help='gas target (default for gaseous elements and gas compounds)')
    g.add_argument('--solid', dest='gas', action='store_false', help='solid target')
    p.add_argument('--mass', type=float, metavar='AMU', help='ion mass (default: most abundant isotope)')
    p.add_argument('--correction', type=float, default=None,
                   help='compound (Bragg) correction factor (default 1)')


# ----------------------------------------------------------------- commands
def cmd_run(a):
    from .trimin import TrimInput
    from .srim import run_trim
    from . import sr
    if len(a.files) > 1 and a.output is not None:
        raise UserError('-o/--output can be used with one input file only')
    for path in a.files:
        if not os.path.isfile(path):
            raise UserError('no such file: %s' % path)
        with open(path, 'rb') as fh:
            head = fh.read(40).decode('latin-1')
        if head.startswith('---Stopping'):
            out = sr.run(path, a.output)
            if not a.quiet:
                print('stopping table written to %s' % out)
            continue
        if len(a.files) > 1:
            outdir = os.path.join(os.path.dirname(os.path.abspath(path)), 'SRIM Outputs')
        else:
            outdir = 'SRIM Outputs' if a.output is None else a.output
        ti = TrimInput.read(path)
        n = int(ti.number_ions if a.number_ions is None else a.number_ions)
        t0 = time.time()
        state = {'last': 0.0}

        def progress(i, eng):
            now = time.time()
            if now - state['last'] > 1 or i == n:
                state['last'] = now
                sys.stderr.write('\rion %d / %d  (%.1f s)' % (i, n, now - t0))
                sys.stderr.flush()

        eng = run_trim(ti, outdir, number_ions=n, progress=None if a.quiet else progress)
        if not a.quiet:
            tal = eng.tal
            sys.stderr.write('\n')
            print('ions: %d  backscattered: %d  transmitted: %d' % (eng.ion_number, tal.n_back, tal.n_trans))
            print('average range %.6g A, straggling %.6g A' % (tal.g_a38, tal.g_a48))
            print('output written to %s' % os.path.abspath(outdir))
    return 0


def cmd_stopping(a):
    kw = dict(emin=api.parse_energy(a.emin), emax=api.parse_energy(a.emax), units=api.stopping_units(a.units),
              density=a.density, gas=a.gas, mass=a.mass, correction=a.correction)
    if a.srim:
        t = api.stopping_report(a.ion, a.target, **kw)
        if a.output:
            t.write(a.output)
        else:
            sys.stdout.write('\n'.join(t.lines()) + '\n')
        return 0
    z, m = api.resolve_ion(a.ion, a.mass)
    mat = api.material(a.target, a.density, a.gas, a.correction)
    rows = api.stopping_table(a.ion, a.target, **kw)
    units = api.STOPPING_UNITS[kw['units']]
    if a.csv:
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator='\n')
        w.writerow(['energy_keV', 'S_elec [%s]' % units, 'S_nuc [%s]' % units, 'range_A',
                    'long_straggling_A', 'lat_straggling_A'])
        for r in rows:
            w.writerow(['%.6g' % r[0]] + ['%.6e' % v for v in r[1:3]] + ['%.6g' % v for v in r[3:]])
        _write(buf.getvalue(), a.output)
        return 0
    o = ['Ion:    %s' % _ion_text(z, m),
         'Target: %s' % _material_text(mat),
         'Stopping units: %s' % units,
         '',
         '%-12s %11s %11s %12s %12s %12s' % ('Energy', 'dE/dx Elec', 'dE/dx Nucl', 'Projected', 'Longitudinal',
                                               'Lateral'),
         '%-12s %11s %11s %12s %12s %12s' % ('', '', '', 'Range', 'Straggling', 'Straggling'),
         '-' * 75]
    for e, se, sn, rp, lon, lat in rows:
        o.append('%-12s %11.3E %11.3E %12s %12s %12s' % (_fmt_energy(e), se, sn, _fmt_length(rp),
                                                          _fmt_length(lon), _fmt_length(lat)))
    _write('\n'.join(o) + '\n', a.output)
    return 0


def cmd_eloss(a):
    value, unit = api.parse_thickness(a.thickness)
    energy = api.parse_energy(a.energy)
    r = api.energy_loss(a.ion, a.target, energy=energy, thickness=value, unit=unit, angle=a.angle,
                        density=a.density, gas=a.gas, mass=a.mass, correction=a.correction)
    mat = r.pop('material')
    if a.json:
        r['target'] = mat.as_dict()
        r['ion'] = dict(zip(('z', 'mass'), api.resolve_ion(a.ion, a.mass)))
        if not a.profile:
            r.pop('profile')
        print(json.dumps(r, indent=1))
        return 0
    z, m = api.resolve_ion(a.ion, a.mass)
    o = ['Ion:    %s, %s' % (_ion_text(z, m), _fmt_energy(energy)),
         'Target: %s' % _material_text(mat),
         'Layer:  %s = %s' % (a.thickness, _fmt_length(r['thickness'] / r['path_factor'])) +
         ('' if not a.angle else ', at %g deg: path %s' % (a.angle, _fmt_length(r['thickness']))),
         '        %.5g atoms/cm2, %.5g ug/cm2%s' % (r['areal_density'], r['mass_thickness'],
                                                 ' along the path' if a.angle else ''),
         '']
    if r['stopped']:
        o.append('The ions stop in the layer (mean range about %s).' % _fmt_length(r['depth_reached']))
    else:
        o += ['exit energy        %10.6g keV' % r['energy_out'],
              'energy loss        %10.6g keV  (electronic %.6g, nuclear %.6g)'
              % (r['energy_loss'], r['electronic_loss'], r['nuclear_loss']),
              'mean energy        %10.6g keV  (average over the layer)' % r['mean_energy'],
              'straggling sigma   %10.5g keV  FWHM %.5g keV  (TRIM model)' % (r['sigma'], r['fwhm']),
              'straggling sigma   %10.5g keV  FWHM %.5g keV  (Bohr)' % (r['sigma_bohr'], r['fwhm_bohr'])]
    print('\n'.join(o))
    return 0


def cmd_trim(a):
    value, unit = api.parse_thickness(a.thickness)
    energy = api.parse_energy(a.energy)
    z, m = api.resolve_ion(a.ion, a.mass)
    mat = api.material(a.target, a.density, a.gas, a.correction)
    n = a.ions
    t0 = time.time()
    state = {'last': 0.0}

    def progress(i, eng):
        now = time.time()
        if now - state['last'] > 1 or i == n:
            state['last'] = now
            sys.stderr.write('\rion %d / %d  (%.1f s)' % (i, n, now - t0))
            sys.stderr.flush()

    res = api.trim(a.ion, mat, energy=energy, thickness=value, unit=unit, ions=n, calculation=a.calculation,
                   angle=a.angle, seed=a.seed, directory=a.out, mass=a.mass,
                   progress=None if a.quiet else progress)
    if not a.quiet:
        sys.stderr.write('\n')
    rg = res.range
    o = ['Ion:    %s, %s, %d ions' % (_ion_text(z, m), _fmt_energy(energy), n),
         'Target: %s, %s' % (_material_text(mat), a.thickness),
         '']
    if rg is not None:
        o.append('mean projected range   %s' % _fmt_length(rg.average_range))
        o.append('straggling             %s' % _fmt_length(rg.straggling))
        if rg.lateral_range is not None:
            o.append('lateral range          %s' % _fmt_length(rg.lateral_range))
    o.append('backscattered ions     %d' % (len(res.backscat) if res.backscat is not None else 0))
    o.append('transmitted ions       %d' % (len(res.transmit) if res.transmit is not None else 0))
    if res.vacancy is not None and res.vacancy.vacancies_per_ion is not None:
        o.append('vacancies per ion      %g' % res.vacancy.vacancies_per_ion)
    o.append('')
    o.append('output files in %s' % res.directory)
    print('\n'.join(o))
    return 0


def cmd_web(a):
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401
    except ImportError:
        raise UserError('the web application needs fastapi and uvicorn: pip install "lunasrim[web]"')
    from .web.__main__ import main as web_main
    return web_main(['--host', a.host, '--port', str(a.port)]) or 0


# ----------------------------------------------------------------- parser
def build_parser():
    p = argparse.ArgumentParser(
        prog='lunasrim',
        description='LunaSRIM %s: SRIM-2013 (TRIM Monte Carlo and stopping/range tables) in Python. '
                    '"lunasrim FILE" runs a TRIM.IN or SR.IN file, as "lunasrim run FILE".' % __version__,
        epilog='Energies: keV, or a number with eV/keV/MeV/GeV (400, 400keV, 2MeV). Thicknesses: a number with '
               'A, nm, um, mm, cm, ug/cm2, mg/cm2, 1e15/cm2 or atoms/cm2 (1um, 200nm, 2e19atoms/cm2). '
               'Run "lunasrim COMMAND -h" for the options of a command.')
    p.add_argument('-V', '--version', action='version', version='lunasrim %s' % __version__)
    sub = p.add_subparsers(dest='command', metavar='COMMAND')

    q = sub.add_parser('run', help='run TRIM.IN or SR.IN files',
                       description='Run SRIM input files: TRIM.IN (TRIM calculation, output files in '
                                   '"SRIM Outputs") or SR.IN (stopping/range table, written next to SR.IN as '
                                   'named in it); the type is detected from the file.')
    q.add_argument('files', nargs='+', metavar='FILE', help='TRIM.IN or SR.IN file')
    q.add_argument('-o', '--output', default=None,
                   help='directory for the output files (default: "SRIM Outputs" for TRIM.IN, the directory of '
                        'SR.IN for SR.IN; with several files, next to each file)')
    q.add_argument('-n', '--number-ions', type=int, default=None,
                   help='number of ions (default: the number given in TRIM.IN)')
    q.add_argument('-q', '--quiet', action='store_true', help='no progress output')
    q.set_defaults(func=cmd_run)

    q = sub.add_parser('stopping', help='stopping and range table',
                       description='Stopping powers, projected range and straggling of an ion in a target '
                                   '(SRIM "Stopping and Range Tables", identical to SRModule.exe).')
    _target_args(q)
    q.add_argument('--emin', default='10keV', help='lowest energy (default 10keV)')
    q.add_argument('--emax', default='10MeV', help='highest energy (default 10MeV)')
    q.add_argument('--units', default='MeV/(mg/cm2)',
                   help='stopping units: 1-8 or eV/A, keV/um, MeV/mm, keV/(ug/cm2), MeV/(mg/cm2), keV/(mg/cm2), '
                        '"eV/(1E15 atoms/cm2)", LSS (default MeV/(mg/cm2))')
    f = q.add_mutually_exclusive_group()
    f.add_argument('--csv', action='store_true', help='CSV output (energy keV, stopping, ranges in A)')
    f.add_argument('--srim', action='store_true', help='SRModule.exe output file format')
    q.add_argument('-o', '--output', metavar='FILE', help='write to FILE instead of standard output')
    q.set_defaults(func=cmd_stopping)

    q = sub.add_parser('eloss', help='energy loss and energy straggling in a layer',
                       description='Mean energy loss and energy straggling (TRIM model and Bohr) of an ion beam '
                                   'crossing one layer, without Monte Carlo.')
    _target_args(q)
    q.add_argument('--energy', required=True, help='ion energy (keV, or e.g. 400keV, 2MeV)')
    q.add_argument('--thickness', required=True,
                   help='layer thickness with unit: 1um, 200nm, 50A, 100ug/cm2, 0.5mg/cm2, 2e19atoms/cm2')
    q.add_argument('--angle', type=float, default=0.0, metavar='DEG',
                   help='beam tilt from the layer normal in degrees (default 0)')
    q.add_argument('--json', action='store_true', help='all results as JSON')
    q.add_argument('--profile', action='store_true', help='with --json: include the depth profile')
    q.set_defaults(func=cmd_eloss)

    q = sub.add_parser('trim', help='TRIM Monte Carlo in one layer',
                       description='TRIM Monte Carlo calculation of an ion beam in one layer, with the TRIM '
                                   'default target energies; writes TRIM.IN and "SRIM Outputs" into --out.')
    _target_args(q)
    q.add_argument('--energy', required=True, help='ion energy (keV, or e.g. 50keV, 2MeV)')
    q.add_argument('--thickness', required=True, help='layer thickness with unit: 200nm, 1um, 100ug/cm2, ...')
    q.add_argument('--ions', type=int, default=1000, help='number of ions (default 1000)')
    q.add_argument('--calculation', type=int, choices=(1, 2, 3), default=1,
                   help='1 = ion distribution and quick damage (Kinchin-Pease, default), '
                        '2 = full damage cascades, 3 = monolayer collision steps / sputtering')
    q.add_argument('--angle', type=float, default=0.0, metavar='DEG', help='angle of incidence (default 0)')
    q.add_argument('--seed', type=int, default=None, help='random number seed (default: random)')
    q.add_argument('--out', default='lunasrim_output', metavar='DIR',
                   help='output directory (default ./lunasrim_output)')
    q.add_argument('-q', '--quiet', action='store_true', help='no progress output')
    q.set_defaults(func=cmd_trim)

    q = sub.add_parser('web', help='start the web application',
                       description='Start the LunaSRIM web application (needs: pip install "lunasrim[web]").')
    q.add_argument('--host', default='127.0.0.1', help='default 127.0.0.1')
    q.add_argument('--port', type=int, default=8000, help='default 8000')
    q.set_defaults(func=cmd_web)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] not in COMMANDS and argv[0] not in ('-h', '--help', '-V', '--version'):
        argv.insert(0, 'run')           # backwards compatible: lunasrim TRIM.IN [-o DIR] [-n N]
    p = build_parser()
    a = p.parse_args(argv)
    if a.command is None:
        p.print_help(sys.stderr)
        return 2
    try:
        return a.func(a)
    except (UserError, ValueError, KeyError, OSError) as ex:
        msg = ex.args[0] if isinstance(ex, KeyError) and ex.args else str(ex)
        sys.stderr.write('lunasrim: error: %s\n' % msg)
        return 2
    except KeyboardInterrupt:
        sys.stderr.write('\ninterrupted\n')
        return 130


if __name__ == '__main__':
    sys.exit(main())
