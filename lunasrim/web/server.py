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
import os
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from typing import List, Optional

from .. import __version__
from . import core

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
MAX_IONS = int(os.environ.get('LUNASRIM_MAX_IONS', '100000'))
MAX_JOBS = int(os.environ.get('LUNASRIM_MAX_JOBS', '2'))
KEEP_JOBS = 40

app = FastAPI(title='LunaSRIM', version=__version__)


def _call(fn, *args):
    try:
        return fn(*args)
    except core.ApiError as ex:
        raise HTTPException(ex.status, ex.detail)


def _dict(q):
    return q.model_dump() if hasattr(q, 'model_dump') else q.dict()


# ----------------------------------------------------------------- data
@app.get('/api/elements')
def api_elements():
    return core.elements()


@app.get('/api/compounds')
def api_compounds():
    return core.compound_list()


@app.get('/api/info')
def api_info():
    return core.info(MAX_IONS)


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


class ElossIn(BaseModel):
    ion: IonIn
    target: MaterialIn
    thickness: float = Field(..., gt=0)
    unit: str = 'nm'               # A, nm, um, mm, ug/cm2, mg/cm2, 1e15/cm2, atoms/cm2


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


# ------------------------------------------------------------ stopping
@app.post('/api/stopping')
def api_stopping(q: StoppingIn):
    return _call(core.stopping, _dict(q))


@app.post('/api/eloss')
def api_eloss(q: ElossIn):
    return _call(core.eloss, _dict(q))


@app.post('/api/implant')
def api_implant(q: ImplantIn):
    return _call(core.implant, _dict(q))


# ----------------------------------------------------------------- TRIM
JOBS = {}
POOL = ThreadPoolExecutor(max_workers=MAX_JOBS)
LOCK = threading.Lock()


def _prune():
    with LOCK:
        jobs = sorted(JOBS.items(), key=lambda kv: kv[1].created)
        while len(jobs) > KEEP_JOBS:
            jid, j = jobs.pop(0)
            if j.status in ('finished', 'cancelled', 'error'):
                j.remove()
                JOBS.pop(jid, None)
            else:
                break


@app.post('/api/trim')
def api_trim(q: TrimIn):
    req = _call(core.trim_request, _dict(q), MAX_IONS)
    _prune()
    job = core.TrimRun(req)
    jid = uuid.uuid4().hex[:12]
    JOBS[jid] = job
    POOL.submit(job.step)
    return {'id': jid}


def _job(jid):
    job = JOBS.get(jid)
    if job is None:
        raise HTTPException(404, 'unknown job')
    return job


@app.get('/api/trim/{jid}')
def api_trim_status(jid: str, live: int = 0):
    out = _job(jid).state(bool(live))
    out['id'] = jid
    return out


@app.post('/api/trim/{jid}/cancel')
def api_trim_cancel(jid: str):
    job = _job(jid)
    job.cancel = True
    return {'id': jid, 'status': job.status}


@app.get('/api/trim/{jid}/files.zip')
def api_trim_zip(jid: str):
    data = _call(_job(jid).zip_bytes)
    return Response(data, media_type='application/zip',
                    headers={'Content-Disposition': 'attachment; filename="lunasrim_%s.zip"' % jid})


@app.get('/api/trim/{jid}/file/{name}')
def api_trim_file(jid: str, name: str):
    p = _call(_job(jid).file_path, name)
    return FileResponse(p, media_type='text/plain', filename=name)


# --------------------------------------------------------------- static
@app.get('/vendor/plotly.min.js')
def plotly_js():
    p = os.path.join(STATIC, 'vendor', 'plotly.min.js')
    if os.path.exists(p):
        return FileResponse(p, media_type='application/javascript')
    return Response(status_code=404)


app.mount('/', StaticFiles(directory=STATIC, html=True), name='static')
