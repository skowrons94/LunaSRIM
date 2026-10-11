// LunaSRIM in the browser: Pyodide runs the lunasrim package in this worker.
// Messages: {id, method, path, body} -> {id, status, text} (JSON text, as from
// the HTTP server), {id, kind: 'zip'|'file', job, name} -> {id, status, bytes|text}.
let py = null, api = null, ready = null;
const busy = new Set();

function init(base) {
  if (!ready) ready = (async () => {
    importScripts(base + 'pyodide/pyodide.js');
    py = await loadPyodide({ indexURL: base + 'pyodide/' });
    const r = await fetch(base + 'lunasrim.zip');
    if (!r.ok) throw new Error('cannot load lunasrim.zip (' + r.status + ')');
    py.unpackArchive(await r.arrayBuffer(), 'zip', { extractDir: '/home/pyodide/lib' });
    py.runPython("import sys; sys.path.insert(0, '/home/pyodide/lib')");
    api = py.pyimport('lunasrim.web.pyodide_api');
  })();
  return ready;
}

const pause = () => new Promise(r => setTimeout(r, 0));

async function runTrim(job) {
  // the calculation runs in slices so that progress and cancel requests are answered
  if (busy.has(job)) return;
  busy.add(job);
  try {
    while (api.trim_step(job, 0.2)) await pause();
  } finally { busy.delete(job); }
}

self.onmessage = async (ev) => {
  const m = ev.data;
  try {
    await init(m.base);
    if (m.kind === 'zip') {
      const b = api.trim_zip(m.job).toJs();
      self.postMessage({ id: m.id, status: 200, bytes: b }, [b.buffer]);
      return;
    }
    if (m.kind === 'file') {
      self.postMessage({ id: m.id, status: 200, text: api.trim_file(m.job, m.name) });
      return;
    }
    const res = api.handle(m.method, m.path, m.body || '');
    const [status, text] = res.toJs();
    res.destroy();
    self.postMessage({ id: m.id, status, text });
    if (m.method === 'POST' && m.path === '/api/trim' && status === 200) runTrim(JSON.parse(text).id);
  } catch (e) {
    self.postMessage({ id: m.id, status: 500, text: JSON.stringify({ detail: String(e && e.message || e) }) });
  }
};
