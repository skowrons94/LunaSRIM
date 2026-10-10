"""LunaSRIM web application: HTTP API on top of the lunasrim package.

Run with ``python -m lunasrim.web`` (needs ``fastapi`` and ``uvicorn``).

Endpoints (JSON):

* ``GET  /api/elements``          element data (masses, densities, default energies)
* ``GET  /api/compounds``         SRIM's compound dictionary
* ``POST /api/stopping``          stopping / range / straggling table (SR module)
* ``POST /api/implant``           projected range, straggling and implantation profiles
* ``POST /api/trim``              start a TRIM Monte Carlo calculation (returns a job id)
* ``GET  /api/trim/{id}``         progress and, when finished, the results
* ``POST /api/trim/{id}/cancel``  stop a calculation (results of the ions done so far)
* ``GET  /api/trim/{id}/files.zip``  TRIM.IN and all output files
"""
import io
import math
import os
import shutil
import tempfile
import threading
import time
import traceback
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from typing import List, Optional

from .. import __version__, compounds, data_files, sr
from ..output import TrimOutput, TrimDistributions
from ..results import Results
from ..srim import Element, Ion, Layer, Target, TRIM
from ..stopping import GAS_ELEMENTS
from ..trim import TrimEngine
from ..trimin import TrimInput

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
MAX_IONS = int(os.environ.get('LUNASRIM_MAX_IONS', '100000'))
MAX_JOBS = int(os.environ.get('LUNASRIM_MAX_JOBS', '2'))
KEEP_JOBS = 40

# default displacement / lattice energies proposed by the TRIM user interface
DISPLACEMENT = {1: 10, 6: 28, 7: 28, 8: 28, 14: 15}

app = FastAPI(title='LunaSRIM', version=__version__)


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
    raise HTTPException(400, 'unknown element Z=%s' % z)


@app.get('/api/elements')
def api_elements():
    return elements()


@app.get('/api/compounds')
def api_compounds():
    return [c.as_dict() for c in compounds.load()]


@app.get('/api/info')
def api_info():
    return {'version': __version__, 'max_ions': MAX_IONS, 'units': sr.UNITS}


# ------------------------------------------------------------- models
class IonIn(BaseModel):
    z: int = Field(..., ge=1, le=92)
    mass: Optional[float] = None
    energy: float = 100.0          # keV
    angle: float = 0.0             # degrees


class AtomIn(BaseModel):
    z: int = Field(..., ge=1, le=92)
    stoich: float = 1.0
    mass: Optional[float] = None
    displacement: Optional[float] = None
    lattice: Optional[float] = None
    surface: Optional[float] = None


class MaterialIn(BaseModel):
    name: str = 'target'
    elements: List[AtomIn]
    density: float = Field(..., gt=0)
    gas: bool = False
    correction: float = 1.0        # compound (Bragg) correction
    width: float = 10000.0         # Angstrom (TRIM layers)


class StoppingIn(BaseModel):
    ion: IonIn
    target: MaterialIn
    emin: float = 10.0             # keV
    emax: float = 10000.0
    units: int = Field(5, ge=1, le=8)


class ImplantIn(BaseModel):
    ion: IonIn
    target: MaterialIn
    energies: List[float]          # keV
    doses: Optional[List[float]] = None   # ions/cm2 per energy
    depth_points: int = 400


class TrimIn(BaseModel):
    ion: IonIn
    layers: List[MaterialIn]
    number_ions: int = Field(1000, ge=1)
    calculation: int = Field(1, ge=1, le=3)
    seed: Optional[int] = None
    trajectories: int = Field(50, ge=0, le=500)
    plot_xmin: float = 0.0        # Angstrom, depth window of the distributions
    plot_xmax: float = 0.0        # 0 = whole target
    exyz: float = Field(0.0, ge=0)          # EXYZ.txt energy interval (eV), 0 = off
    collisions: int = Field(0, ge=0, le=2)  # COLLISON.txt: 1 = ion collisions, 2 = with cascades
    title: str = 'LunaSRIM web'


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
        raise HTTPException(400, 'the target has no elements')
    si.elements = els
    si.units = units
    si.emin, si.emax = emin, emax
    return si


def _table(si):
    t = sr.StoppingTable(si, strict=False)
    try:
        t.calculate()
    except sr.SrError as ex:
        raise HTTPException(422, str(ex))
    except (ValueError, ZeroDivisionError, OverflowError) as ex:
        raise HTTPException(422, 'calculation failed: %s' % ex)
    return t


# ------------------------------------------------------------ stopping
@app.post('/api/stopping')
def api_stopping(q: StoppingIn):
    if not (0 < q.emin <= q.emax):
        raise HTTPException(400, 'need 0 < E-min <= E-max')
    if q.emax / max(q.emin, 1e-12) > 1e12 or q.emax / _ion_mass(q.ion) > 1e7:
        raise HTTPException(400, 'energy range too large (max 12 decades, 10 GeV/u)')
    si = _sr_input(q.ion, q.target, q.units, q.emin, q.emax)
    t = _table(si)
    u = t.unit_factor[q.units] / t.conv
    rows = [{'energy': e, 'se': u * se, 'sn': u * sn, 'range': rp, 'long': lon, 'lat': lat}
            for (e, se, sn, rp, lon, lat) in t.rows]
    text = '\r\n'.join(t.lines()) + '\r\n'
    factors = {k: t.unit_factor[k] / t.unit_factor[q.units] for k in range(1, 9)}
    return {'rows': rows, 'units': sr.UNITS[q.units].strip(), 'unit_factors': factors,
            'atom_density': t.atom_density, 'sr_output': text,
            'composition': [{'z': e.z, 'fraction': float(f)} for e, f in zip(si.elements, t.frac)]}


# ------------------------------------------------------------ implantation
@app.post('/api/implant')
def api_implant(q: ImplantIn):
    if not q.energies:
        raise HTTPException(400, 'no energies')
    if len(q.energies) > 50:
        raise HTTPException(400, 'at most 50 energies')
    res = []
    for k, en in enumerate(q.energies):
        if en <= 0:
            raise HTTPException(400, 'energies must be positive')
        si = _sr_input(q.ion, q.target, 5, en, en)
        t = _table(si)
        e, se, sn, rp, lon, lat = t.rows[-1]
        u = t.unit_factor[5] / t.conv
        dose = q.doses[k] if q.doses and k < len(q.doses) else 1e15
        res.append({'energy': en, 'se': u * se, 'sn': u * sn, 'range': rp,
                    'long': lon, 'lat': lat, 'dose': dose})
    # Gaussian depth profiles (atoms/cm3), depth in Angstrom
    xmax = max(r['range'] + 4 * r['long'] for r in res)
    npt = max(50, min(q.depth_points, 2000))
    x = [xmax * i / (npt - 1) for i in range(npt)]
    comps = []
    total = [0.0] * npt
    for r in res:
        s = max(r['long'], 1e-6)
        c = [r['dose'] / (math.sqrt(2 * math.pi) * s * 1e-8) *
             math.exp(-0.5 * ((xi - r['range']) / s) ** 2) for xi in x]
        comps.append(c)
        total = [a + b for a, b in zip(total, c)]
    n_target = (6.022e23 * q.target.density) / _mean_mass(q.target)
    return {'results': res, 'depth': x, 'profiles': comps, 'total': total,
            'target_atom_density': n_target}


def _mean_mass(m):
    tot = sum(a.stoich for a in m.elements) or 1.0
    return sum(a.stoich / tot * (a.mass or element(a.z)['mass']) for a in m.elements)


# ----------------------------------------------------------------- TRIM
class Job(object):
    def __init__(self, q):
        self.id = uuid.uuid4().hex[:12]
        self.q = q
        self.dir = tempfile.mkdtemp(prefix='lunasrim_job_')
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


JOBS = {}
POOL = ThreadPoolExecutor(max_workers=MAX_JOBS)
LOCK = threading.Lock()


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
            raise HTTPException(400, 'layer %s has no elements' % m.name)
        if not (m.width > 0):
            raise HTTPException(400, 'layer widths must be positive')
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


def _run_job(job):
    try:
        job.status = 'running'
        job.started = time.time()
        trim = _trim_input(job.q)
        text = trim.trim_input_text()
        with open(os.path.join(job.dir, 'TRIM.IN'), 'wb') as fh:
            fh.write(text.encode('latin-1'))
        ti = TrimInput.parse(text)
        eng = TrimEngine(ti)
        eng.out = TrimOutput(eng, job.dir)
        eng.out.start()
        ntraj = job.q.trajectories
        try:
            for i in range(job.total):
                if job.cancel:
                    break
                eng.record_path = i < ntraj
                st = eng.run_ion()
                if eng.record_path:
                    pts = _thin([(round(x, 3), round(y, 3), round(z, 3)) for (x, y, z) in eng.path])
                    job.trajectories.append({'fate': eng.fate, 'points': pts})
                    if eng.recoil_segments and len(job.recoils) < 60000:
                        job.recoils.extend([tuple(round(v, 3) for v in s) for s in eng.recoil_segments[:3000]])
                if eng.fate == 'S' and len(job.live['range']) < 20000:
                    job.live['range'].append(st.x)
                job.live['back'] = eng.tal.n_back
                job.live['trans'] = eng.tal.n_trans
                job.done = i + 1
        finally:
            eng.out.close()
        if job.done == 0:
            raise RuntimeError('no ion was calculated')
        TrimDistributions(eng, job.dir).write_all()
        job.result = _collect(job, eng)
        job.status = 'cancelled' if job.cancel else 'finished'
    except HTTPException as ex:
        job.status, job.error = 'error', ex.detail
    except Exception as ex:          # pragma: no cover - reported to the client
        job.status = 'error'
        job.error = '%s: %s' % (type(ex).__name__, ex)
        traceback.print_exc()
    finally:
        job.finished = time.time()


def _lst(a):
    return [float(v) for v in a] if a is not None else []


def _collect(job, eng):
    r = Results(job.dir)
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
        out['ioniz'] = {'depth': _lst(r.ioniz.depth), 'ions': _lst(r.ioniz.ions), 'recoils': _lst(r.ioniz.recoils)}
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
    out['trajectories'] = job.trajectories
    out['recoil_segments'] = job.recoils
    out['files'] = sorted(os.listdir(job.dir))
    return out


def _prune():
    with LOCK:
        jobs = sorted(JOBS.values(), key=lambda j: j.created)
        while len(jobs) > KEEP_JOBS:
            j = jobs.pop(0)
            if j.status in ('finished', 'cancelled', 'error'):
                shutil.rmtree(j.dir, ignore_errors=True)
                JOBS.pop(j.id, None)
            else:
                break


@app.post('/api/trim')
def api_trim(q: TrimIn):
    if q.number_ions > MAX_IONS:
        raise HTTPException(400, 'at most %d ions on this server' % MAX_IONS)
    if not q.layers:
        raise HTTPException(400, 'no target layers')
    if len(q.layers) > 100:
        raise HTTPException(400, 'at most 100 layers')
    if not (0 < q.ion.energy <= 1e7):
        raise HTTPException(400, 'ion energy must be between 0 and 10 GeV')
    if not (-89.99 <= q.ion.angle <= 89.99):
        raise HTTPException(400, 'angle must be within +-89.99 degrees')
    if q.plot_xmax and not (0 <= q.plot_xmin < q.plot_xmax):
        raise HTTPException(400, 'the depth window needs 0 <= from < to')
    _trim_input(q)                      # validate before queueing
    _prune()
    job = Job(q)
    JOBS[job.id] = job
    POOL.submit(_run_job, job)
    return {'id': job.id}


def _job(jid):
    job = JOBS.get(jid)
    if job is None:
        raise HTTPException(404, 'unknown job')
    return job


@app.get('/api/trim/{jid}')
def api_trim_status(jid: str, live: int = 0):
    job = _job(jid)
    now = time.time()
    el = ((job.finished or now) - job.started) if job.started else 0.0
    out = {'id': job.id, 'status': job.status, 'done': job.done, 'total': job.total,
           'elapsed': el, 'error': job.error}
    if job.status in ('finished', 'cancelled'):
        out['result'] = job.result
    elif live:
        rng = job.live['range']
        out['live'] = {'range': rng[-5000:], 'backscattered': job.live['back'],
                       'transmitted': job.live['trans'],
                       'trajectories': job.trajectories[-20:]}
    return out


@app.post('/api/trim/{jid}/cancel')
def api_trim_cancel(jid: str):
    job = _job(jid)
    job.cancel = True
    return {'id': job.id, 'status': job.status}


@app.get('/api/trim/{jid}/files.zip')
def api_trim_zip(jid: str):
    job = _job(jid)
    if job.status not in ('finished', 'cancelled'):
        raise HTTPException(409, 'the calculation is not finished')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
        for name in sorted(os.listdir(job.dir)):
            z.write(os.path.join(job.dir, name), name)
    buf.seek(0)
    return StreamingResponse(buf, media_type='application/zip',
                             headers={'Content-Disposition': 'attachment; filename="lunasrim_%s.zip"' % jid})


@app.get('/api/trim/{jid}/file/{name}')
def api_trim_file(jid: str, name: str):
    job = _job(jid)
    if '/' in name or '\\' in name or name.startswith('.'):
        raise HTTPException(400, 'bad file name')
    p = os.path.join(job.dir, name)
    if not os.path.isfile(p):
        raise HTTPException(404, 'no such file')
    return FileResponse(p, media_type='text/plain', filename=name)


# --------------------------------------------------------------- static
@app.get('/vendor/plotly.min.js')
def plotly_js():
    p = os.path.join(STATIC, 'vendor', 'plotly.min.js')
    if os.path.exists(p):
        return FileResponse(p, media_type='application/javascript')
    return Response(status_code=404)


app.mount('/', StaticFiles(directory=STATIC, html=True), name='static')
