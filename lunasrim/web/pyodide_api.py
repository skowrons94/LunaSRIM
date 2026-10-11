"""The web application's API inside the browser (Pyodide).

``static/pyworker.js`` loads the package into Pyodide and calls
:func:`handle` with the same method/path/JSON body that the FastAPI server
receives; the answer is ``(status, json_text)``.  A TRIM calculation is run
in slices by :func:`trim_step`, so that the worker can answer progress and
cancel requests in between.
"""
import json
import re

from . import core

MAX_IONS = 100000
RUNS = {}
_next = [0]


def _ok(data):
    return 200, json.dumps(data)


def handle(method, path, body_text):
    try:
        body = json.loads(body_text) if body_text else {}
        path, _, query = path.partition('?')
        if method == 'GET' and path == '/api/elements':
            return _ok(core.elements())
        if method == 'GET' and path == '/api/compounds':
            return _ok(core.compound_list())
        if method == 'GET' and path == '/api/info':
            d = core.info(MAX_IONS)
            d['static'] = True
            return _ok(d)
        if method == 'POST' and path == '/api/stopping':
            return _ok(core.stopping(body))
        if method == 'POST' and path == '/api/implant':
            return _ok(core.implant(body))
        if method == 'POST' and path == '/api/trim':
            q = core.trim_request(body, MAX_IONS)
            for jid in list(RUNS):            # keep the memory file system small
                if RUNS[jid].status in ('finished', 'cancelled', 'error'):
                    RUNS.pop(jid).remove()
            _next[0] += 1
            jid = 'run%d' % _next[0]
            RUNS[jid] = core.TrimRun(q)
            return _ok({'id': jid})
        m = re.match(r'^/api/trim/([\w-]+)(/cancel)?$', path)
        if m:
            run = RUNS.get(m.group(1))
            if run is None:
                raise core.ApiError(404, 'unknown job')
            if m.group(2):
                run.cancel = True
                return _ok({'id': m.group(1), 'status': run.status})
            out = run.state('live=1' in query)
            out['id'] = m.group(1)
            return _ok(out)
        raise core.ApiError(404, 'not found: %s %s' % (method, path))
    except core.ApiError as ex:
        return ex.status, json.dumps({'detail': ex.detail})
    except Exception as ex:
        return 500, json.dumps({'detail': '%s: %s' % (type(ex).__name__, ex)})


def trim_step(jid, seconds):
    """advance a TRIM calculation; True while it is not over"""
    run = RUNS.get(jid)
    return bool(run is not None and run.step(seconds))


def trim_zip(jid):
    return RUNS[jid].zip_bytes()


def trim_file(jid, name):
    return RUNS[jid].file_text(name)
