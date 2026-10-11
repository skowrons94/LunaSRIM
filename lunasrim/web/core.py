"""Calculations behind the LunaSRIM web application, without any web framework.

The same functions serve the FastAPI server (:mod:`lunasrim.web.server`) and
the static site, where they run in the browser with Pyodide
(``static/pyworker.js``).  Requests and answers are plain JSON-like dicts;
errors are raised as :class:`ApiError` (an HTTP status and a message).
"""
import io
import math
import os
import shutil
import tempfile
import time
import zipfile
from types import SimpleNamespace

from .. import __version__, compounds, data_files, sr, straggling
from ..output import TrimOutput, TrimDistributions
from ..results import Results
from ..srim import Element, Ion, Layer, Target, TRIM
from ..stopping import GAS_ELEMENTS
from ..trim import TrimEngine
from ..trimin import TrimInput

# default displacement / lattice energies proposed by the TRIM user interface
DISPLACEMENT = {1: 10, 6: 28, 7: 28, 8: 28, 14: 15}


class ApiError(Exception):
    def __init__(self, status, detail):
        Exception.__init__(self, detail)
        self.status = status
        self.detail = detail


# ------------------------------------------------------------------- input
_LABEL = {'z': 'element', 'emin': 'minimum energy', 'emax': 'maximum energy', 'energy': 'ion energy',
          'stoich': 'stoichiometry', 'density': 'density', 'width': 'layer width',
          'number_ions': 'number of ions'}


def _num(d, key, default=None, kind=float, lo=None, hi=None, lo_open=False, required=False):
    """a number of the request: ``default`` if the key is missing; null (an
    empty or invalid input field in the browser) gives the default too,
    or an error if the value is ``required``"""
    v = d.get(key, default) if isinstance(d, dict) else default
    if v is None:
        if required:
            raise ApiError(422, '%s: enter a number' % _LABEL.get(key, key))
        return default
    try:
        if kind is int:
            if isinstance(v, float) and v != int(v):
                raise ValueError
            v = int(v)
        else:
            v = float(v)
    except (TypeError, ValueError):
        raise ApiError(422, '%s: not a valid number' % key)
    if not math.isfinite(v):
        raise ApiError(422, '%s: not a valid number' % key)
    if lo is not None and (v <= lo if lo_open else v < lo):
        raise ApiError(422, '%s: must be %s %s' % (key, '>' if lo_open else '>=', lo))
    if hi is not None and v > hi:
        raise ApiError(422, '%s: must be <= %s' % (key, hi))
    return v


def _ion(d):
    if not isinstance(d, dict):
        raise ApiError(422, 'ion: field required')
    return SimpleNamespace(z=_num(d, 'z', kind=int, lo=1, hi=92, required=True), mass=_num(d, 'mass'),
                           energy=_num(d, 'energy'), angle=_num(d, 'angle', 0.0))


def _atom(d):
    if not isinstance(d, dict):
        raise ApiError(422, 'elements: not valid')
    return SimpleNamespace(z=_num(d, 'z', kind=int, lo=1, hi=92, required=True),
                           stoich=_num(d, 'stoich', 1.0, required=True),
                           mass=_num(d, 'mass'), displacement=_num(d, 'displacement'),
                           lattice=_num(d, 'lattice'), surface=_num(d, 'surface'))


def _material(d):
    if not isinstance(d, dict):
        raise ApiError(422, 'target: field required')
    els = d.get('elements')
    if not isinstance(els, list):
        raise ApiError(422, 'elements: field required')
    return SimpleNamespace(name=str(d.get('name', 'target') or ''), elements=[_atom(a) for a in els],
                           density=_num(d, 'density', lo=0, lo_open=True, required=True),
                           gas=bool(d.get('gas', False)), correction=_num(d, 'correction', 1.0),
                           width=_num(d, 'width', 10000.0, required=True))


# ----------------------------------------------------------------- data
def _elements():
    atoms = data_files.load_atomdata()
    scoef = data_files.load_scoef()
    out = []
    for z in sorted(atoms):
        a = atoms[z]
        sc = scoef.get(z)
        gas = z in GAS_ELEMENTS
        out.append({
            'z': z, 'symbol': a.symbol, 'name': a.name,
            'mass': a.natural_weight, 'ion_mass': a.mai_weight,
            'density': (sc.gas_density if (gas and sc) else a.density),
            'solid_density': a.density,
            'gas': gas,
            'displacement': DISPLACEMENT.get(z, 25), 'lattice': 3.0,
            'surface': a.heat_sublimation,
        })
    return out


_ELEMENTS = None


def elements():
    global _ELEMENTS
    if _ELEMENTS is None:
        _ELEMENTS = _elements()
    return _ELEMENTS


def element(z):
    for e in elements():
        if e['z'] == z:
            return e
    raise ApiError(400, 'unknown element Z=%s' % z)


def compound_list():
    return [c.as_dict() for c in compounds.load()]


def info(max_ions):
    return {'version': __version__, 'max_ions': max_ions, 'units': sr.UNITS}


# ------------------------------------------------------------ stopping
def _ion_mass(ion):
    return ion.mass if ion.mass else element(ion.z)['ion_mass']


def _sr_input(ion, target, units, emin, emax):
    si = sr.SrInput()
    si.output_name = ('%s in %s' % (element(ion.z)['name'], target.name or 'target'))[:60]
    si.ion_z = ion.z
    si.ion_mass = _ion_mass(ion)
    si.phase = 1 if target.gas else 0
    si.density = target.density
    si.correction = target.correction
    els = []
    for a in target.elements:
        if a.stoich <= 0:
            continue
        e = element(a.z)
        els.append(sr.SrElement(a.z, e['name'], a.stoich, a.mass if a.mass else e['mass']))
    if not els:
        raise ApiError(400, 'the target has no elements')
    si.elements = els
    si.units = units
    si.emin, si.emax = emin, emax
    return si


def _table(si):
    t = sr.StoppingTable(si, strict=False)
    try:
        t.calculate()
    except sr.SrError as ex:
        raise ApiError(422, str(ex))
    except (ValueError, ZeroDivisionError, OverflowError) as ex:
        raise ApiError(422, 'calculation failed: %s' % ex)
    return t


def stopping(d):
    """POST /api/stopping"""
    ion, target = _ion(d.get('ion')), _material(d.get('target'))
    emin, emax = _num(d, 'emin', 10.0, required=True), _num(d, 'emax', 10000.0, required=True)
    units = _num(d, 'units', 5, kind=int, lo=1, hi=8)
    if not (0 < emin <= emax):
        raise ApiError(400, 'need 0 < E-min <= E-max')
    if emax / max(emin, 1e-12) > 1e12 or emax / _ion_mass(ion) > 1e7:
        raise ApiError(400, 'energy range too large (max 12 decades, 10 GeV/u)')
    si = _sr_input(ion, target, units, emin, emax)
    t = _table(si)
    u = t.unit_factor[units] / t.conv
    rows = [{'energy': e, 'se': u * se, 'sn': u * sn, 'range': rp, 'long': lon, 'lat': lat}
            for (e, se, sn, rp, lon, lat) in t.rows]
    text = '\r\n'.join(t.lines()) + '\r\n'
    factors = {k: t.unit_factor[k] / t.unit_factor[units] for k in range(1, 9)}
    return {'rows': rows, 'units': sr.UNITS[units].strip(), 'unit_factors': factors,
            'atom_density': t.atom_density, 'sr_output': text,
            'composition': [{'z': e.z, 'fraction': float(f)} for e, f in zip(si.elements, t.frac)]}


# ------------------------------------------------------------ implantation
def _mean_mass(m):
    tot = sum(a.stoich for a in m.elements) or 1.0
    return sum(a.stoich / tot * (a.mass or element(a.z)['mass']) for a in m.elements)


def implant(d):
    """POST /api/implant"""
    ion, target = _ion(d.get('ion')), _material(d.get('target'))
    energies = d.get('energies')
    doses = d.get('doses')
    npt = _num(d, 'depth_points', 400, kind=int)
    if not isinstance(energies, list) or not energies:
        raise ApiError(400, 'no energies')
    if len(energies) > 50:
        raise ApiError(400, 'at most 50 energies')
    energies = [_num({'energy': v}, 'energy', required=True) for v in energies]
    res = []
    for k, en in enumerate(energies):
        if en <= 0:
            raise ApiError(400, 'energies must be positive')
        si = _sr_input(ion, target, 5, en, en)
        t = _table(si)
        e, se, sn, rp, lon, lat = t.rows[-1]
        u = t.unit_factor[5] / t.conv
        dose = float(doses[k]) if doses and k < len(doses) and doses[k] is not None else 1e15
        res.append({'energy': en, 'se': u * se, 'sn': u * sn, 'range': rp,
                    'long': lon, 'lat': lat, 'dose': dose})
    # Gaussian depth profiles (atoms/cm3), depth in Angstrom
    xmax = max(r['range'] + 4 * r['long'] for r in res)
    npt = max(50, min(npt, 2000))
    x = [xmax * i / (npt - 1) for i in range(npt)]
    comps = []
    total = [0.0] * npt
    for r in res:
        s = max(r['long'], 1e-6)
        c = [r['dose'] / (math.sqrt(2 * math.pi) * s * 1e-8) *
             math.exp(-0.5 * ((xi - r['range']) / s) ** 2) for xi in x]
        comps.append(c)
        total = [a + b for a, b in zip(total, c)]
    n_target = (6.022e23 * target.density) / _mean_mass(target)
    return {'results': res, 'depth': x, 'profiles': comps, 'total': total,
            'target_atom_density': n_target}


# -------------------------------------------------- energy loss / straggling
def eloss(d):
    """POST /api/eloss: energy loss and straggling of ions crossing a layer"""
    ion, target = _ion(d.get('ion')), _material(d.get('target'))
    if ion.energy is None:
        raise ApiError(422, 'ion energy: enter a number')
    thickness = _num(d, 'thickness', required=True, lo=0, lo_open=True)
    unit = str(d.get('unit') or 'nm')
    if unit not in straggling.THICKNESS_UNITS:
        raise ApiError(422, 'unit: one of %s' % ', '.join(straggling.THICKNESS_UNITS))
    if not (0 < ion.energy <= 1e7):
        raise ApiError(400, 'ion energy must be between 0 and 10 GeV')
    els = []
    for a in target.elements:
        if a.stoich > 0:
            els.append((a.z, a.stoich, a.mass if a.mass else element(a.z)['mass']))
    if not els:
        raise ApiError(400, 'the target has no elements')
    if not (-89.99 <= ion.angle <= 89.99):
        raise ApiError(400, 'angle must be within +-89.99 degrees')
    tilt = 1 / math.cos(math.radians(ion.angle))     # path length / layer thickness
    try:
        r = straggling.energy_loss(ion.z, _ion_mass(ion), ion.energy, els, target.density, thickness * tilt, unit,
                                   phase=1 if target.gas else 0, correction=target.correction)
    except straggling.StragglingError as ex:
        raise ApiError(400, str(ex))
    except (ValueError, ZeroDivisionError, OverflowError) as ex:
        raise ApiError(422, 'calculation failed: %s' % ex)
    prof = r.pop('profile')
    n = len(prof['depth'])
    step = max(1, n // 300)
    idx = list(range(0, n, step))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    r['profile'] = {k: [v[i] for i in idx] for k, v in prof.items()}
    r['angle'], r['path_factor'] = ion.angle, tilt
    return r


# ----------------------------------------------------------------- TRIM
def trim_request(d, max_ions):
    """validated TRIM request (POST /api/trim)"""
    q = SimpleNamespace(
        ion=_ion(d.get('ion')),
        layers=[_material(m) for m in (d.get('layers') or [])],
        number_ions=_num(d, 'number_ions', 1000, kind=int, lo=1, required=True),
        calculation=_num(d, 'calculation', 1, kind=int, lo=1, hi=3),
        seed=_num(d, 'seed', kind=int),
        trajectories=_num(d, 'trajectories', 50, kind=int, lo=0, hi=500),
        plot_xmin=_num(d, 'plot_xmin', 0.0), plot_xmax=_num(d, 'plot_xmax', 0.0),
        exyz=_num(d, 'exyz', 0.0, lo=0), collisions=_num(d, 'collisions', 0, kind=int, lo=0, hi=2),
        title=str(d.get('title', 'LunaSRIM web') or ''))
    if q.ion.energy is None:
        raise ApiError(422, 'ion energy: enter a number')
    if q.number_ions > max_ions:
        raise ApiError(400, 'at most %d ions here' % max_ions)
    if not q.layers:
        raise ApiError(400, 'no target layers')
    if len(q.layers) > 100:
        raise ApiError(400, 'at most 100 layers')
    if not (0 < q.ion.energy <= 1e7):
        raise ApiError(400, 'ion energy must be between 0 and 10 GeV')
    if not (-89.99 <= q.ion.angle <= 89.99):
        raise ApiError(400, 'angle must be within +-89.99 degrees')
    if q.plot_xmax and not (0 <= q.plot_xmin < q.plot_xmax):
        raise ApiError(400, 'the depth window needs 0 <= from < to')
    _trim_input(q)                      # validate the target
    return q


def _trim_input(q):
    ion = Ion(q.ion.z, energy=q.ion.energy * 1000.0, mass=_ion_mass(q.ion))
    layers = []
    for m in q.layers:
        els = {}
        for a in m.elements:
            if a.stoich <= 0:
                continue
            e = element(a.z)
            key = Element(a.z, a.mass if a.mass else e['mass'])
            els[key] = {'stoich': a.stoich,
                        'E_d': e['displacement'] if a.displacement is None else a.displacement,
                        'lattice': e['lattice'] if a.lattice is None else a.lattice,
                        'surface': e['surface'] if a.surface is None else a.surface}
        if not els:
            raise ApiError(400, 'layer %s has no elements' % m.name)
        if not (m.width > 0):
            raise ApiError(400, 'layer widths must be positive')
        layers.append(Layer(els, density=m.density, width=m.width, phase=1 if m.gas else 0,
                            name=(m.name or 'layer')[:20].replace('"', ''),
                            bragg_correction=m.correction))
    seed = q.seed if q.seed is not None else int(time.time()) % 100000
    trim = TRIM(Target(layers), ion, number_ions=q.number_ions, calculation=q.calculation,
                angle_ions=q.ion.angle, random_seed=seed, ranges=True, backscattered=True,
                plot_xmin=q.plot_xmin if q.plot_xmax > 0 else 0, plot_xmax=q.plot_xmax,
                transmit=True, sputtered=True, collisions=q.collisions, exyz=q.exyz,
                description=q.title[:60] or 'LunaSRIM web')
    return trim


def _thin(points, nmax=400):
    if len(points) <= nmax:
        return points
    step = len(points) / float(nmax - 1)
    out = [points[int(i * step)] for i in range(nmax - 1)]
    out.append(points[-1])
    return out


def _thin_curve(c, nmax=400):
    if len(c['x']) > nmax:
        idx = _thin(list(range(len(c['x']))), nmax)
        c['x'] = [c['x'][i] for i in idx]
        c['e'] = [c['e'][i] for i in idx]
    return c


def _lst(a):
    return [float(v) for v in a] if a is not None else []


class TrimRun(object):
    """One TRIM calculation, advanced ion by ion with :meth:`step` (the
    server runs it in a thread, the browser in slices between messages)."""

    def __init__(self, q, directory=None):
        self.q = q
        self.dir = directory or tempfile.mkdtemp(prefix='lunasrim_job_')
        self.status = 'queued'
        self.done = 0
        self.total = q.number_ions
        self.started = None
        self.finished = None
        self.cancel = False
        self.error = None
        self.result = None
        self.trajectories = []
        self.recoils = []
        self.live = {'range': [], 'back': 0, 'trans': 0}
        self.created = time.time()
        self.eng = None

    def _start(self):
        self.status = 'running'
        self.started = time.time()
        trim = _trim_input(self.q)
        text = trim.trim_input_text()
        with open(os.path.join(self.dir, 'TRIM.IN'), 'wb') as fh:
            fh.write(text.encode('latin-1'))
        eng = TrimEngine(TrimInput.parse(text))
        eng.out = TrimOutput(eng, self.dir)
        eng.out.start()
        self.eng = eng

    def step(self, seconds=None):
        """run ions for about ``seconds`` (all of them if None); returns
        True while the calculation is not over"""
        if self.status in ('finished', 'cancelled', 'error'):
            return False
        try:
            if self.eng is None:
                self._start()
            eng = self.eng
            t_end = None if seconds is None else time.time() + seconds
            ntraj = self.q.trajectories
            while self.done < self.total and not self.cancel:
                i = self.done
                eng.record_path = i < ntraj
                st = eng.run_ion()
                if eng.record_path:
                    pts = _thin([(round(x, 3), round(y, 3), round(z, 3)) for (x, y, z) in eng.path])
                    self.trajectories.append({'fate': eng.fate, 'points': pts})
                    if eng.recoil_segments and len(self.recoils) < 60000:
                        self.recoils.extend([tuple(round(v, 3) for v in s) for s in eng.recoil_segments[:3000]])
                if eng.fate == 'S' and len(self.live['range']) < 20000:
                    self.live['range'].append(st.x)
                self.live['back'] = eng.tal.n_back
                self.live['trans'] = eng.tal.n_trans
                self.done = i + 1
                if t_end is not None and time.time() > t_end:
                    return True
            self._finish()
        except ApiError as ex:
            self._fail(ex.detail)
        except Exception as ex:          # reported to the client
            self._fail('%s: %s' % (type(ex).__name__, ex))
        return False

    def _fail(self, msg):
        if self.eng is not None and self.eng.out is not None:
            self.eng.out.close()
        self.status, self.error = 'error', msg
        self.finished = time.time()

    def _finish(self):
        eng = self.eng
        eng.out.close()
        if self.done == 0:
            raise RuntimeError('no ion was calculated')
        TrimDistributions(eng, self.dir).write_all()
        self.result = self._collect()
        self.status = 'cancelled' if self.cancel else 'finished'
        self.finished = time.time()

    def state(self, live=False):
        """GET /api/trim/{id}"""
        now = time.time()
        el = ((self.finished or now) - self.started) if self.started else 0.0
        out = {'status': self.status, 'done': self.done, 'total': self.total, 'elapsed': el,
               'error': self.error}
        if self.status in ('finished', 'cancelled'):
            out['result'] = self.result
        elif live:
            out['live'] = {'range': self.live['range'][-5000:], 'backscattered': self.live['back'],
                           'transmitted': self.live['trans'], 'trajectories': self.trajectories[-20:]}
        return out

    def zip_bytes(self):
        if self.status not in ('finished', 'cancelled'):
            raise ApiError(409, 'the calculation is not finished')
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
            for name in sorted(os.listdir(self.dir)):
                z.write(os.path.join(self.dir, name), name)
        return buf.getvalue()

    def file_path(self, name):
        if '/' in name or '\\' in name or name.startswith('.'):
            raise ApiError(400, 'bad file name')
        p = os.path.join(self.dir, name)
        if not os.path.isfile(p):
            raise ApiError(404, 'no such file')
        return p

    def file_text(self, name):
        with open(self.file_path(name), 'rb') as fh:
            return fh.read().decode('latin-1')

    def remove(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _collect(self):
        eng = self.eng
        r = Results(self.dir)
        t = eng.t
        out = {'ions': eng.ion_number,
               'layers': [{'name': l.name, 'top': (t.depth[i - 1] if i else 0.0), 'bottom': t.depth[i]}
                          for i, l in enumerate(eng.ti.layers)],
               'elements': [{'z': e.z, 'symbol': e.symbol} for e in eng.ti.elements],
               'summary': {}}
        s = out['summary']
        if r.range is not None:
            s.update(range=r.range.average_range, straggling=r.range.straggling,
                     lateral=r.range.lateral_range, lateral_straggling=r.range.lateral_straggling,
                     radial=r.range.radial_range, radial_straggling=r.range.radial_straggling)
            out['range'] = {'depth': _lst(r.range.depth), 'ions': _lst(r.range.ions),
                            'recoils': [_lst(c) for c in r.range.elements]}
        if r.vacancy is not None:
            s.update(vacancies=r.vacancy.vacancies_per_ion, displacements=r.vacancy.displacements_per_ion,
                     replacements=r.vacancy.replacements_per_ion)
            out['vacancy'] = {'depth': _lst(r.vacancy.depth), 'knock_ons': _lst(r.vacancy.knock_ons),
                              'vacancies': [_lst(c) for c in r.vacancy.vacancies]}
        s.update(backscattered=eng.tal.n_back, transmitted=eng.tal.n_trans,
                 sputtered=len(r.sputter) if r.sputter is not None else 0)
        if r.ioniz is not None:
            out['ioniz'] = {'depth': _lst(r.ioniz.depth), 'ions': _lst(r.ioniz.ions),
                            'recoils': _lst(r.ioniz.recoils)}
        if r.phonons is not None:
            out['phonons'] = {'depth': _lst(r.phonons.depth), 'ions': _lst(r.phonons.ions),
                              'recoils': _lst(r.phonons.recoils)}
        if r.e2recoil is not None:
            out['e2recoil'] = {'depth': _lst(r.e2recoil.depth), 'ions': _lst(r.e2recoil.ions),
                               'absorbed': [_lst(c) for c in r.e2recoil.absorbed]}
        if r.lateral is not None:
            out['lateral'] = {'depth': _lst(r.lateral.depth), 'lateral_range': _lst(r.lateral.lateral_range),
                              'lateral_straggling': _lst(r.lateral.lateral_straggling),
                              'radial_range': _lst(r.lateral.radial_range),
                              'radial_straggling': _lst(r.lateral.radial_straggling)}
        if r.novac is not None:
            out['novac'] = {'depth': _lst(r.novac.depth), 'number': _lst(r.novac.number)}

        def atoms(f, nmax=20000):
            if f is None or len(f) == 0:
                return None
            n = min(len(f), nmax)
            return {'energy': _lst(f.energy[:n]), 'x': _lst(f.position[0][:n]), 'y': _lst(f.position[1][:n]),
                    'z': _lst(f.position[2][:n]), 'cx': _lst(f.direction[0][:n]),
                    'cy': _lst(f.direction[1][:n]), 'cz': _lst(f.direction[2][:n]),
                    'atom': _lst(f.z[:n]), 'count': len(f)}
        out['transmitted'] = atoms(r.transmit)
        out['backscattered'] = atoms(r.backscat)
        out['sputtered'] = atoms(r.sputter)
        if r.range3d is not None and len(r.range3d):
            n = min(len(r.range3d), 20000)
            out['stopped'] = {'x': _lst(r.range3d.positions[0][:n]), 'y': _lst(r.range3d.positions[1][:n]),
                              'z': _lst(r.range3d.positions[2][:n]), 'count': len(r.range3d)}
        if r.exyz is not None and len(r.exyz):
            curves = {}
            for n, e, x in zip(r.exyz.ion_number, r.exyz.energy, r.exyz.positions[0]):
                n = int(n)
                if n not in curves:
                    if len(curves) >= 100:
                        break
                    curves[n] = {'n': n, 'x': [], 'e': []}
                curves[n]['x'].append(float(x))
                curves[n]['e'].append(float(e))
            out['exyz'] = [_thin_curve(c) for c in curves.values()]
        out['trajectories'] = self.trajectories
        out['recoil_segments'] = self.recoils
        out['files'] = sorted(os.listdir(self.dir))
        return out
