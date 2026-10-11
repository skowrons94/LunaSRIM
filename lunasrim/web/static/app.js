/* LunaSRIM web application */
'use strict';

const S = { elements: [], byZ: {}, bySym: {}, compounds: [], info: {} };
const $ = (s, r) => (r || document).querySelector(s);
const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (v == null || v === false) continue;
    if (k === 'class') e.className = v;
    else if (k === 'html') e.innerHTML = v;
    else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v === true ? '' : v);
  }
  for (const k of kids.flat()) if (k != null && k !== false) e.append(k.nodeType ? k : document.createTextNode(k));
  return e;
}

function toast(msg, err) {
  const t = $('#toast');
  t.textContent = msg;
  t.className = 'toast' + (err ? ' err' : '');
  clearTimeout(toast._t);
  toast._t = setTimeout(() => t.classList.add('hidden'), err ? 6000 : 2500);
}

// Static site (GitHub Pages): no server, lunasrim runs in the browser with
// Pyodide in a web worker (pyworker.js) that answers the same /api/ requests.
const LS = window.LUNASRIM_STATIC ? {
  w: null, n: 0, pend: {},
  base: new URL('.', location.href).href,
  send(msg) {
    if (!this.w) {
      this.w = new Worker('pyworker.js');
      this.w.onmessage = (ev) => { const p = this.pend[ev.data.id]; delete this.pend[ev.data.id]; if (p) p(ev.data); };
    }
    const id = ++this.n;
    return new Promise(res => { this.pend[id] = res; this.w.postMessage(Object.assign({ id, base: this.base }, msg)); });
  },
} : null;

async function api(path, body) {
  if (LS) {
    const r = await LS.send({ method: body === undefined ? 'GET' : 'POST', path, body: body === undefined ? '' : JSON.stringify(body) });
    const data = JSON.parse(r.text);
    if (r.status !== 200) throw new Error(data && data.detail ? String(data.detail) : 'error ' + r.status);
    return data;
  }
  const opt = body === undefined ? {} : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
  const r = await fetch(path, opt);
  let data = null;
  try { data = await r.json(); } catch (e) { /* not json */ }
  if (!r.ok) {
    let msg = r.statusText;
    if (data && data.detail) msg = typeof data.detail === 'string' ? data.detail : data.detail.map(d => (d.loc || []).slice(1).join('.') + ': ' + d.msg).join('; ');
    throw new Error(msg);
  }
  return data;
}

// ------------------------------------------------------------- formatting
function fmt(v, d = 4) {
  if (v == null || !isFinite(v)) return '–';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(d - 1).replace('e+', 'E+').replace('e-', 'E-');
  return Number(v.toPrecision(d)).toString();
}
function fmtLenAs(A, ref, d = 4) {  // in the unit chosen for ref
  if (A == null || !isFinite(A)) return '–';
  const [u, f] = lenUnit(Math.abs(ref || A) * 2);
  return fmt(A / f, d) + ' ' + u;
}
function fmtLen(A, d = 4) {      // Angstrom -> best unit
  if (A == null || !isFinite(A)) return '–';
  const a = Math.abs(A);
  if (a < 1e3) return fmt(A, d) + ' Å';
  if (a < 1e4) return fmt(A / 10, d) + ' nm';
  if (a < 1e7) return fmt(A / 1e4, d) + ' µm';
  if (a < 1e10) return fmt(A / 1e7, d) + ' mm';
  return fmt(A / 1e10, d) + ' m';
}
function fmtE(keV) {
  if (keV < 1) return fmt(keV * 1000, 4) + ' eV';
  if (keV < 1000) return fmt(keV, 4) + ' keV';
  if (keV < 1e6) return fmt(keV / 1000, 4) + ' MeV';
  return fmt(keV / 1e6, 4) + ' GeV';
}
const EUNIT = { eV: 1e-3, keV: 1, MeV: 1e3, GeV: 1e6 };
function lenUnit(maxA, pref) {
  if (pref && pref !== 'auto') return { A: ['Å', 1], nm: ['nm', 10], um: ['µm', 1e4] }[pref];
  if (maxA < 2e3) return ['Å', 1];
  if (maxA < 2e4) return ['nm', 10];
  if (maxA < 2e7) return ['µm', 1e4];
  return ['mm', 1e7];
}

// ------------------------------------------------------------------ plots
// Plots in the style of a journal figure: closed frame, ticks inside (minor
// ticks too), faint grid, thin lines.  Series colours are a fixed categorical
// order checked for colour-vision deficiencies (light and dark steps in
// style.css, --s1..--s8); lines of one plot also differ in dash pattern.
const C = (i) => cssVar('--s' + (i % 8 + 1));
const DASHES = ['solid', 'dash', 'dot', 'dashdot', 'longdash', 'longdashdot', '6px,2px,1px,2px', '1px,3px'];
function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
function rgba(hex, a) {
  const n = parseInt(hex.replace('#', ''), 16);
  return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`;
}
const isNarrow = () => window.innerWidth < 640;
const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif';

function plot(id, traces, layout = {}, opts = {}) {
  const el = document.getElementById(id);
  if (!window.Plotly) { el.textContent = 'Plotting library not loaded (no network?). Run tools/fetch_web_vendor.sh on the server.'; return; }
  el._plot = [traces, layout, opts];
  el.classList.remove('cleared');
  // a hidden plot is drawn when its tab is shown (resizePlots): drawn now it
  // would get Plotly's default width and overflow a phone screen
  if (!el.offsetParent) { el._stale = true; return; }
  el._stale = false;
  const narrow = isNarrow();
  const text = cssVar('--text'), muted = cssVar('--muted'), grid = cssVar('--grid'), panel = cssVar('--panel'), line = cssVar('--line');
  const fs = narrow ? 11 : 12.5;
  // decades of positive data on an axis (labels only at powers of ten when wide)
  const decades = (key) => {
    let lo = Infinity, hi = 0;
    for (const t of traces) for (const v of (t[key] || [])) if (v > 0 && isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    return hi > 0 ? Math.log10(hi / lo) : 0;
  };
  const ax = (a, key) => {
    a = Object.assign({}, a || {});
    const title = typeof a.title === 'string' ? { text: a.title } : (a.title || {});
    delete a.title;
    const log = a.type === 'log';
    const wide = log && (a.range ? a.range[1] - a.range[0] : decades(key)) >= 2;
    return Object.assign({
      showline: true, mirror: 'ticks', linecolor: text, linewidth: 1,
      ticks: 'inside', ticklen: 6, tickwidth: 1, tickcolor: text,
      minor: { ticks: 'inside', ticklen: 3, tickwidth: 1, tickcolor: text, showgrid: false, dtick: wide ? 'D1' : undefined },
      showgrid: true, gridcolor: grid, gridwidth: 1, zeroline: false,
      tickfont: { color: text, size: fs }, automargin: true,
      exponentformat: 'power', dtick: wide ? 1 : undefined,
      title: Object.assign({ font: { size: fs + 0.5, color: text }, standoff: 10 }, title),
    }, a);
  };
  // thin lines; dash patterns as a second cue when a plot has several lines
  let k = 0;
  const lines = traces.filter(t => (t.mode || '').includes('lines') && !(t.line && t.line.width === 0) && t.showlegend !== false).length;
  traces = traces.map(t => {
    t = Object.assign({}, t);
    if ((t.mode || '').includes('lines')) {
      t.line = Object.assign({ width: 1.8 }, t.line);
      if (lines > 1 && t.line.width && t.showlegend !== false && !t.line.dash) t.line.dash = DASHES[k++ % DASHES.length];
    }
    return t;
  });
  const L = Object.assign({
    paper_bgcolor: panel, plot_bgcolor: panel, font: { family: FONT, color: text, size: fs },
    margin: narrow ? { l: 8, r: 10, t: 8, b: 8 } : { l: 12, r: 18, t: 12, b: 12 },
    legend: { orientation: 'h', x: 0, xanchor: 'left', y: 1.01, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)', font: { size: fs }, itemwidth: 40, tracegroupgap: 4 },
    hovermode: 'x unified',
    hoverlabel: { bgcolor: panel, bordercolor: line, font: { family: FONT, color: text, size: fs } },
    dragmode: narrow ? false : 'zoom',
  }, layout);
  L.xaxis = ax(layout.xaxis, 'x'); L.yaxis = ax(layout.yaxis, 'y');
  // large or small values on a linear y axis: plot y / 10^n, "(×10^n)" in the title
  if (L.yaxis.type !== 'log' && !layout.noScale) {
    let m = 0;
    for (const t of traces) if (t.y && t.type !== 'histogram2d') for (const v of t.y) if (isFinite(v)) m = Math.max(m, Math.abs(v));
    const n = m > 0 ? Math.floor(Math.log10(m)) : 0;
    if (n >= 4 || n <= -3) {
      const f = Math.pow(10, -n);
      traces = traces.map(t => (t.y && t.type !== 'histogram2d') ? Object.assign({}, t, { y: t.y.map(v => v * f) }) : t);
      const tt = L.yaxis.title.text || '';
      L.yaxis.title = Object.assign({}, L.yaxis.title, { text: `${tt}${tt ? ' ' : ''}(×10<sup>${n}</sup>)` });
      if (L.yaxis.range) L.yaxis.range = L.yaxis.range.map(v => v * f);
    }
  }
  if (layout.title) {
    const t = typeof layout.title === 'string' ? layout.title : layout.title.text;
    L.title = { text: t, x: 0, xref: 'paper', xanchor: 'left', y: 1, yref: 'container', yanchor: 'top', pad: { t: 6 }, font: { size: fs + 1.5, color: text }, automargin: true };
  }
  if (L.shapes) L.shapes = L.shapes.map(s => Object.assign({ layer: 'below' }, s));
  Plotly.react(el, traces, L, Object.assign({
    responsive: true, displaylogo: false, scrollZoom: false,
    displayModeBar: narrow ? false : 'hover',
    modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d', 'zoomIn2d', 'zoomOut2d', 'toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'],
    toImageButtonOptions: { format: 'svg', filename: id },
  }, opts));
}
function clearPlot(id, html) {
  const el = document.getElementById(id);
  if (window.Plotly) Plotly.purge(el);
  el._plot = null; el.innerHTML = html; el.classList.add('cleared');
}
function resizePlots(root) {
  if (!window.Plotly) return;
  // a plot drawn while hidden has no size: lay it out again
  $$('.plot', root).forEach(p => { if (p._plot && p.offsetParent) plot(p.id, ...p._plot); });
}
// phone <-> desktop layout: draw the plots again
let _narrow = null;
window.addEventListener('resize', () => {
  const n = isNarrow();
  if (_narrow !== null && n !== _narrow) $$('.plot').forEach(p => { if (p._plot) plot(p.id, ...p._plot); });
  _narrow = n;
});
// histogram of a list of values as a step line (bin centres)
function histogram(values, lo, hi, n) {
  const w = (hi - lo) / n || 1, c = new Array(n).fill(0);
  for (const v of values) { const i = Math.min(n - 1, Math.max(0, Math.floor((v - lo) / w))); c[i]++; }
  return { x: c.map((_, i) => lo + (i + 0.5) * w), y: c };
}
const stepTrace = (x, y, name, color, extra) => Object.assign({ x, y, name, mode: 'lines', line: { shape: 'hvh', color, width: 1.5 }, fill: 'tozeroy', fillcolor: rgba(color, 0.12) }, extra || {});

// ----------------------------------------------------- periodic table
const POS = {};        // Z -> [row, col]
(function () {
  const rows = [[1, 2], [3, 10], [11, 18], [19, 36], [37, 54], [55, 86], [87, 118]];
  const put = (z, r, c) => { POS[z] = [r, c]; };
  put(1, 1, 1); put(2, 1, 18);
  for (const [r, a] of [[2, 3], [3, 11]]) { put(a, r, 1); put(a + 1, r, 2); for (let i = 0; i < 6; i++) put(a + 2 + i, r, 13 + i); }
  for (const [r, a] of [[4, 19], [5, 37]]) for (let i = 0; i < 18; i++) put(a + i, r, 1 + i);
  for (const [r, a, lan] of [[6, 55, 57], [7, 87, 89]]) {
    put(a, r, 1); put(a + 1, r, 2);
    for (let i = 0; i < 15; i++) put(lan + i, r + 3, 3 + i);      // lanthanides / actinides
    for (let i = 0; i < 15; i++) put(lan + 15 + i, r, 4 + i);
  }
  void rows;
})();
function category(z) {
  if ([3, 11, 19, 37, 55, 87].includes(z)) return 'alkali';
  if ([4, 12, 20, 38, 56, 88].includes(z)) return 'alkaline';
  if (z >= 57 && z <= 71) return 'lanth';
  if (z >= 89 && z <= 103) return 'actin';
  if ([5, 14, 32, 33, 51, 52].includes(z)) return 'metalloid';
  if ([9, 17, 35, 53, 85].includes(z)) return 'halogen';
  if ([2, 10, 18, 36, 54, 86].includes(z)) return 'noble';
  if ([1, 6, 7, 8, 15, 16, 34].includes(z)) return 'nonmetal';
  if ([13, 31, 49, 50, 81, 82, 83, 84].includes(z)) return 'post';
  return 'transition';
}

const Picker = {
  cb: null,
  build() {
    const g = $('#pt-grid');
    g.innerHTML = '';
    const cells = [];
    for (const e of S.elements) {
      const [r, c] = POS[e.z];
      const cell = h('button', { class: 'pt-cell cat-' + category(e.z), title: `${e.name} (Z=${e.z}, ${e.mass} u)`, style: `grid-row:${r};grid-column:${c}`, 'data-z': e.z,
        onclick: () => { const cb = Picker.cb; $('#pt-dialog').close(); if (cb) cb(e.z); } },
        h('span', { class: 'z' }, e.z), h('span', { class: 's' }, e.symbol), h('span', { class: 'm' }, fmt(e.mass, 4)));
      cells.push(cell);
    }
    g.append(...cells);
    g.append(h('div', { style: 'grid-row:8;grid-column:1/span 18;height:8px' }));
    $('#pt-search').addEventListener('input', (ev) => {
      const q = ev.target.value.trim().toLowerCase();
      $$('.pt-cell', g).forEach(c => {
        const e = S.byZ[+c.dataset.z];
        const hit = !q || e.symbol.toLowerCase() === q || e.name.toLowerCase().startsWith(q) || String(e.z) === q;
        c.classList.toggle('dim', !hit);
      });
    });
    $('#pt-search').addEventListener('keydown', (ev) => {
      if (ev.key !== 'Enter') return;
      const hit = $$('.pt-cell:not(.dim)', g);
      if (hit.length) hit[0].click();
    });
    $('#pt-close').onclick = () => $('#pt-dialog').close();
  },
  open(title, cb) {
    this.cb = cb;
    $('#pt-title').textContent = title;
    $('#pt-search').value = '';
    $$('.pt-cell').forEach(c => c.classList.remove('dim'));
    $('#pt-dialog').showModal();
    setTimeout(() => $('#pt-search').focus(), 30);
  },
};

// ------------------------------------------------------------ ion picker
class IonPicker {
  constructor(root, opts = {}) {
    this.root = root; this.opts = opts;
    this.z = opts.z || 1; this.mass = null; this.energy = opts.energy || 100; this.angle = 0;
    this.render();
  }
  render() {
    const e = S.byZ[this.z];
    this.root.innerHTML = '';
    const tile = h('div', { class: 'el-tile cat-' + category(this.z), title: 'Choose the ion', onclick: () => Picker.open('Select the ion', z => { this.z = z; this.mass = null; this.render(); this.changed(); }) },
      h('span', { class: 'z' }, this.z), h('span', { class: 'sym' }, e.symbol), h('span', { class: 'nm' }, e.name));
    const mass = h('input', { type: 'number', step: 'any', min: 0.5, value: this.mass || e.ion_mass, title: 'Default: mass of the most abundant isotope' });
    mass.addEventListener('change', () => { this.mass = parseFloat(mass.value) || null; this.changed(); });
    const fields = h('div', { class: 'ion-fields' }, h('label', null, 'Mass (amu)', mass));
    if (this.opts.trim) {
      const en = h('input', { type: 'number', step: 'any', min: 0, value: this.energy });
      const eu = h('select', null, ...['eV', 'keV', 'MeV', 'GeV'].map(u => h('option', { selected: u === 'keV' }, u)));
      const set = () => { this.energy = (parseFloat(en.value) || 0) * EUNIT[eu.value]; };
      en.addEventListener('input', set); eu.addEventListener('change', set);
      this._en = en; this._eu = eu;
      const an = h('input', { type: 'number', step: 'any', min: -89.99, max: 89.99, value: this.angle, title: 'Angle of incidence to the surface normal' });
      an.addEventListener('input', () => { this.angle = parseFloat(an.value) || 0; });
      fields.append(h('label', null, 'Energy', h('span', { class: 'unit-input' }, en, eu)), h('label', null, 'Angle (deg)', an));
    } else {
      fields.append(h('label', null, 'Element', h('input', { value: e.name, disabled: true })));
    }
    this.root.append(h('div', { class: 'ion-card' }, tile, fields));
  }
  changed() { if (this.onchange) this.onchange(); }
  get() {
    return { z: this.z, mass: this.mass || S.byZ[this.z].ion_mass, energy: this.energy, angle: this.angle };
  }
  label() { return S.byZ[this.z].symbol; }
}

// ---------------------------------------------------------- target editor
function parseFormula(f) {
  // e.g. SiO2, Al2O3, C8H8, (CH2)n -> {Z: count}
  const out = {};
  const re = /([A-Z][a-z]?)(\d*\.?\d*)/g;
  let m, any = false;
  const s = f.replace(/\s+/g, '');
  while ((m = re.exec(s))) {
    const e = S.bySym[m[1]];
    if (!e) throw new Error('unknown element ' + m[1]);
    out[e.z] = (out[e.z] || 0) + (m[2] ? parseFloat(m[2]) : 1);
    any = true;
  }
  if (!any) throw new Error('empty formula');
  return out;
}

class TargetEditor {
  constructor(root, opts = {}) {
    this.root = root; this.opts = opts;
    this.rows = [];               // {z, stoich, mass, displacement, lattice, surface}
    this.density = null; this.densityAuto = true; this.gas = false; this.correction = 1;
    this.name = opts.name || '';
    this.build();
    this.setElement(opts.z || 14);
  }
  build() {
    const r = this.root;
    r.innerHTML = '';
    this.preset = h('select', { title: 'Element or compound from SRIM\'s compound dictionary' });
    this.preset.append(h('option', { value: '' }, '— element or compound preset —'));
    const cats = {};
    S.compounds.forEach((c, i) => { (cats[c.category] = cats[c.category] || []).push([c, i]); });
    for (const [cat, list] of Object.entries(cats)) {
      const g = h('optgroup', { label: cat });
      list.sort((a, b) => a[0].title.localeCompare(b[0].title)).forEach(([c, i]) => g.append(h('option', { value: i }, c.title)));
      this.preset.append(g);
    }
    this.preset.addEventListener('change', () => { if (this.preset.value !== '') this.setCompound(S.compounds[+this.preset.value]); });
    const addEl = h('button', { class: 'small', title: 'Add an element', onclick: () => Picker.open('Add a target element', z => { this.addRow(z, 1); this.render(); }) }, '+ Element');
    this.advBtn = h('button', { class: 'small', title: 'Masses and (for TRIM) displacement, lattice and surface binding energies', onclick: () => { this.table.classList.toggle('show-adv'); this.advBtn.classList.toggle('primary'); } }, 'Advanced');
    r.append(h('div', { class: 'target-top' }, this.preset));
    const formula = h('input', { placeholder: 'Formula, e.g. SiO2, GaAs, C2H4' });
    const setF = () => {
      try {
        const f = parseFormula(formula.value);
        this.rows = []; for (const [z, n] of Object.entries(f)) this.addRow(+z, n);
        this.name = formula.value.trim(); this.densityAuto = Object.keys(f).length === 1; this.gas = false;
        if (Object.keys(f).length > 1) { this.density = null; toast('Enter the density of ' + this.name); }
        this.preset.value = ''; this.render();
      } catch (e) { toast(e.message, true); }
    };
    formula.addEventListener('keydown', ev => { if (ev.key === 'Enter') setF(); });
    r.append(h('div', { class: 'formula-row' }, formula, h('button', { class: 'small', onclick: setF }, 'Set'), addEl, this.advBtn));
    this.table = h('table', { class: 'comp-table' });
    r.append(this.table);
    this.dens = h('input', { type: 'number', step: 'any', min: 0 });
    this.dens.addEventListener('input', () => { this.density = parseFloat(this.dens.value) || null; this.densityAuto = false; });
    this.gasBox = h('input', { type: 'checkbox' });
    this.gasBox.addEventListener('change', () => { this.gas = this.gasBox.checked; if (this.densityAuto) this.autoDensity(); });
    this.corr = h('input', { type: 'number', step: 'any', min: 0, value: 1, title: 'Compound correction (Bragg rule correction factor), 1 = none' });
    this.corr.addEventListener('input', () => { this.correction = parseFloat(this.corr.value) || 1; });
    r.append(h('div', { class: 'tgt-props' },
      h('label', null, 'Density (g/cm³)', this.dens),
      h('label', null, 'Compound corr.', this.corr),
      h('label', { class: 'chk' }, this.gasBox, 'Gas')));
  }
  addRow(z, stoich) {
    const e = S.byZ[z];
    const ex = this.rows.find(r => r.z === z);
    if (ex) { ex.stoich += stoich; return; }
    this.rows.push({ z, stoich, mass: e.mass, displacement: e.displacement, lattice: e.lattice, surface: e.surface });
  }
  setElement(z) {
    const e = S.byZ[z];
    this.rows = []; this.addRow(z, 1);
    this.name = e.name; this.gas = e.gas; this.densityAuto = true; this.correction = 1;
    this.autoDensity(); this.render();
  }
  setCompound(c) {
    this.rows = [];
    const min = Math.min(...c.elements.map(e => e.fraction));
    for (const e of c.elements) this.addRow(e.z, Number((e.fraction / min).toPrecision(5)));
    this.name = c.name; this.density = c.density; this.densityAuto = false; this.gas = c.gas; this.correction = 1;
    this.render();
  }
  autoDensity() {
    if (this.rows.length === 1) {
      const e = S.byZ[this.rows[0].z];
      this.density = this.gas ? (e.gas ? e.density : e.solid_density * 1e-3) : (e.gas ? e.solid_density : e.density);
    }
  }
  render() {
    const t = this.table;
    const adv = t.classList.contains('show-adv');
    t.innerHTML = '';
    const trim = !!this.opts.trim;
    t.append(h('tr', null, h('th', null, 'Element'), h('th', null, 'Stoich.'), h('th', null, 'Atom %'),
      h('th', { class: 'adv' }, 'Mass'), trim && h('th', { class: 'adv', title: 'Displacement energy (eV)' }, 'E_d'),
      trim && h('th', { class: 'adv', title: 'Lattice binding energy (eV)' }, 'Latt.'), trim && h('th', { class: 'adv', title: 'Surface binding energy (eV)' }, 'Surf.'), h('th')));
    const tot = this.rows.reduce((a, r) => a + r.stoich, 0) || 1;
    this.rows.forEach((row, i) => {
      const e = S.byZ[row.z];
      const num = (k, w) => { const inp = h('input', { type: 'number', step: 'any', min: 0, value: row[k], style: `width:${w}px` }); inp.addEventListener('input', () => { row[k] = parseFloat(inp.value) || 0; if (k === 'stoich') this.updatePct(); }); return inp; };
      const chip = h('span', { class: 'el-chip cat-' + category(row.z), title: e.name + ' — click to change', onclick: () => Picker.open('Replace ' + e.name, z => { const n = S.byZ[z]; Object.assign(row, { z, mass: n.mass, displacement: n.displacement, lattice: n.lattice, surface: n.surface }); if (this.rows.length === 1) { this.name = n.name; this.gas = n.gas; if (this.densityAuto) this.autoDensity(); } this.render(); }) }, e.symbol);
      t.append(h('tr', null, h('td', null, chip), h('td', null, num('stoich', 70)), h('td', { class: 'pct muted' }, fmt(100 * row.stoich / tot, 4)),
        h('td', { class: 'adv' }, num('mass', 74)), trim && h('td', { class: 'adv' }, num('displacement', 52)),
        trim && h('td', { class: 'adv' }, num('lattice', 46)), trim && h('td', { class: 'adv' }, num('surface', 52)),
        h('td', null, this.rows.length > 1 ? h('button', { class: 'x-btn', title: 'Remove', onclick: () => { this.rows.splice(i, 1); if (this.rows.length === 1 && this.densityAuto) this.autoDensity(); this.render(); } }, '✕') : null)));
    });
    this.dens.value = this.density != null ? this.density : '';
    this.gasBox.checked = this.gas;
    this.corr.value = this.correction;
    if (this.onchange) this.onchange();
  }
  updatePct() {
    const tot = this.rows.reduce((a, r) => a + r.stoich, 0) || 1;
    $$('.pct', this.table).forEach((td, i) => { td.textContent = fmt(100 * this.rows[i].stoich / tot, 4); });
  }
  label() {
    if (this.rows.length === 1) return S.byZ[this.rows[0].z].symbol + (this.gas ? ' (gas)' : '');
    return this.name || this.rows.map(r => S.byZ[r.z].symbol + (r.stoich !== 1 ? fmt(r.stoich, 3) : '')).join('');
  }
  get() {
    if (!this.rows.length) throw new Error('add at least one target element');
    if (!(this.density > 0)) throw new Error('enter the target density');
    return {
      name: (this.name || this.label()).slice(0, 40),
      elements: this.rows.map(r => ({ z: r.z, stoich: r.stoich, mass: r.mass, displacement: r.displacement, lattice: r.lattice, surface: r.surface })),
      density: this.density, gas: this.gas, correction: this.correction,
    };
  }
}

// ---------------------------------------------------------------- views
function initTabs() {
  $$('.tab').forEach(b => b.addEventListener('click', () => showView(b.dataset.view)));
  $$('[data-goto]').forEach(a => a.addEventListener('click', ev => { ev.preventDefault(); showView(a.dataset.goto); }));
  $$('.subtabs').forEach(g => $$('.subtab', g).forEach(b => b.addEventListener('click', () => {
    $$('.subtab', g).forEach(x => x.classList.toggle('active', x === b));
    const panel = g.parentElement;
    $$('.pane', panel).forEach(p => p.classList.toggle('active', p.id === b.dataset.pane));
    resizePlots(panel);
    if (b.dataset.pane === 'tr-traj-pane') Trim.draw();
  })));
  const hash = location.hash.slice(1);
  if (hash && $('#view-' + hash)) showView(hash);
}
function showView(v) {
  $$('.tab').forEach(t => t.classList.toggle('active', t.dataset.view === v));
  $$('.view').forEach(x => x.classList.toggle('active', x.id === 'view-' + v));
  history.replaceState(null, '', '#' + v);
  resizePlots($('#view-' + v));
  if (v === 'trim') Trim.draw();
}
function setEmpty(prefix, empty) {
  $('#' + prefix + '-empty').classList.toggle('hidden', !empty);
  $('#' + prefix + '-empty').parentElement.classList.toggle('has-empty', empty);
}
function busy(btn, on, text) {
  if (on) { btn._t = btn.textContent; btn.textContent = text || 'Calculating…'; btn.disabled = true; }
  else { btn.textContent = btn._t || btn.textContent; btn.disabled = false; }
}
function download(name, text, type = 'text/plain') {
  const a = h('a', { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.append(a); a.click(); a.remove();
}

// ------------------------------------------------------- stopping & range
const SR = {
  init() {
    this.ion = new IonPicker($('[data-ion="sr"]'), { z: 1 });
    this.tgt = new TargetEditor($('[data-target="sr"]'), { z: 14 });
    $('#sr-run').onclick = () => this.run();
    $('#sr-csv').onclick = () => this.csv();
    $('#sr-dl').onclick = () => this.data && download((this.title || 'srim') + '.txt', this.data.sr_output);
    setEmpty('sr', true);
  },
  async run() {
    const btn = $('#sr-run');
    let q;
    try {
      q = {
        ion: this.ion.get(), target: this.tgt.get(),
        emin: parseFloat($('#sr-emin').value) * EUNIT[$('#sr-emin-u').value],
        emax: parseFloat($('#sr-emax').value) * EUNIT[$('#sr-emax-u').value],
        units: +$('#sr-units').value,
      };
    } catch (e) { toast(e.message, true); return; }
    busy(btn, true);
    try {
      const d = await api('/api/stopping', q);
      this.data = d; this.q = q;
      this.title = `${this.ion.label()} in ${this.tgt.label()}`;
      setEmpty('sr', false);
      this.render();
    } catch (e) { toast(e.message, true); } finally { busy(btn, false); }
  },
  render() {
    const d = this.data, rows = d.rows, E = rows.map(r => r.energy * 1000);
    const u = d.units;
    plot('plot-sr-stop', [
      { x: E, y: rows.map(r => r.se + r.sn), name: 'Total', mode: 'lines', line: { color: C(2), dash: 'dot', width: 2.2 }, legendrank: 3 },
      { x: E, y: rows.map(r => r.se), name: 'Electronic', mode: 'lines', line: { color: C(0), dash: 'solid' }, legendrank: 1 },
      { x: E, y: rows.map(r => r.sn), name: 'Nuclear', mode: 'lines', line: { color: C(1), dash: 'dash' }, legendrank: 2 },
    ], { title: `Stopping of ${this.title}`, xaxis: { type: 'log', title: 'Ion energy (eV)' }, yaxis: { type: 'log', title: `Stopping power (${u})` } });
    plot('plot-sr-range', [
      { x: E, y: rows.map(r => r.range * 1e-4), name: 'Projected range', mode: 'lines', line: { color: C(0) } },
      { x: E, y: rows.map(r => r.long * 1e-4), name: 'Longitudinal straggling', mode: 'lines', line: { color: C(1) } },
      { x: E, y: rows.map(r => r.lat * 1e-4), name: 'Lateral straggling', mode: 'lines', line: { color: C(2) } },
    ], { title: `Range of ${this.title}`, xaxis: { type: 'log', title: 'Ion energy (eV)' }, yaxis: { type: 'log', title: 'Distance (µm)' } });
    const t = $('#sr-tab');
    t.innerHTML = '';
    t.append(h('tr', null, ...['Energy', `dE/dx elec. (${u})`, `dE/dx nucl. (${u})`, 'Projected range', 'Long. straggling', 'Lat. straggling'].map(x => h('th', null, x))));
    for (const r of rows) t.append(h('tr', null, h('td', null, fmtE(r.energy)), h('td', null, fmt(r.se)), h('td', null, fmt(r.sn)), h('td', null, fmtLen(r.range)), h('td', null, fmtLen(r.long)), h('td', null, fmtLen(r.lat))));
    $('#sr-pre').textContent = d.sr_output;
  },
  csv() {
    if (!this.data) return;
    const u = this.data.units;
    const lines = [`energy_keV,se_${u},sn_${u},range_A,long_straggling_A,lat_straggling_A`];
    for (const r of this.data.rows) lines.push([r.energy, r.se, r.sn, r.range, r.long, r.lat].join(','));
    download((this.title || 'table') + '.csv', lines.join('\n'), 'text/csv');
  },
};

// ------------------------------------------------------------ implantation
const IM = {
  init() {
    this.ion = new IonPicker($('[data-ion="im"]'), { z: 5 });
    this.tgt = new TargetEditor($('[data-target="im"]'), { z: 14 });
    this.rows = [[30, 1e15]];
    this.renderRows();
    $('#im-add').onclick = () => { const last = this.rows[this.rows.length - 1] || [10, 1e15]; this.rows.push([last[0] * 2, last[1]]); this.renderRows(); };
    $('#im-box').onclick = () => this.box();
    $('#im-run').onclick = () => this.run();
    $('#im-log').onchange = () => this.data && this.render();
    $('#im-depthu').onchange = () => this.data && this.render();
    setEmpty('im', true);
  },
  renderRows() {
    const tb = $('#im-energies tbody');
    tb.innerHTML = '';
    this.rows.forEach((r, i) => {
      const e = h('input', { type: 'number', step: 'any', min: 0, value: r[0] });
      const d = h('input', { type: 'text', value: r[1].toExponential(2).replace('e+', 'e') });
      e.oninput = () => { r[0] = parseFloat(e.value) || 0; };
      d.oninput = () => { r[1] = parseFloat(d.value) || 0; };
      tb.append(h('tr', null, h('td', null, e), h('td', null, d), h('td', null, this.rows.length > 1 ? h('button', { class: 'x-btn', onclick: () => { this.rows.splice(i, 1); this.renderRows(); } }, '✕') : null)));
    });
  },
  async box() {
    const spec = prompt('Flat (box) profile: depth from, depth to (in nm) and concentration (atoms/cm³), e.g.\n50, 300, 1e19', '50, 300, 1e19');
    if (!spec) return;
    const [d1, d2, c0] = spec.split(',').map(x => parseFloat(x));
    if (!(d2 > d1 && d1 >= 0 && c0 > 0)) { toast('give: from, to (nm), concentration', true); return; }
    try {
      const t = await api('/api/stopping', { ion: this.ion.get(), target: this.tgt.get(), emin: 0.01, emax: 1e5, units: 5 });
      const rows = t.rows.filter(r => r.range > 0);
      const interp = (x, xs, ys) => { for (let i = 1; i < xs.length; i++) if (xs[i] >= x) { const f = (Math.log(x) - Math.log(xs[i - 1])) / (Math.log(xs[i]) - Math.log(xs[i - 1])); return Math.exp(Math.log(ys[i - 1]) + f * (Math.log(ys[i]) - Math.log(ys[i - 1]))); } return null; };
      const R = rows.map(r => r.range), E = rows.map(r => r.energy), D = rows.map(r => r.long);
      const a = d1 * 10, b = d2 * 10;
      const pts = [];
      let x = Math.max(a, R[0] * 1.01);
      while (x <= b * 1.0001 && pts.length < 20) {
        const s = interp(x, R, D);
        pts.push(x);
        x += Math.max(1.2 * s, (b - a) / 19);
      }
      if (pts.length < 2) pts.push(b);
      const out = [];
      for (let i = 0; i < pts.length; i++) {
        const hstep = i === 0 ? (pts[1] - pts[0]) : i === pts.length - 1 ? (pts[i] - pts[i - 1]) : (pts[i + 1] - pts[i - 1]) / 2;
        const e = interp(pts[i], R, E);
        if (e == null) { toast('depth outside the range table', true); return; }
        out.push([Number(e.toPrecision(4)), Number((c0 * hstep * 1e-8).toPrecision(3))]);
      }
      this.rows = out;
      this.renderRows();
      toast(`${out.length} implants for a flat profile, press Calculate`);
    } catch (e) { toast(e.message, true); }
  },
  async run() {
    const btn = $('#im-run');
    let q;
    try {
      q = { ion: this.ion.get(), target: this.tgt.get(), energies: this.rows.map(r => r[0]), doses: this.rows.map(r => r[1]) };
    } catch (e) { toast(e.message, true); return; }
    busy(btn, true);
    try {
      const [d, c] = await Promise.all([api('/api/implant', q),
        api('/api/stopping', { ion: q.ion, target: q.target, emin: Math.max(0.01, Math.min(...q.energies) / 100), emax: Math.max(...q.energies) * 10, units: 5 })]);
      this.data = d; this.curve = c; this.title = `${this.ion.label()} in ${this.tgt.label()}`;
      setEmpty('im', false);
      this.render();
    } catch (e) { toast(e.message, true); } finally { busy(btn, false); }
  },
  render() {
    const d = this.data;
    const maxd = d.depth[d.depth.length - 1];
    const [lu, ls] = lenUnit(maxd, $('#im-depthu').value);
    const x = d.depth.map(v => v / ls);
    const log = $('#im-log').checked;
    const traces = [];
    if (d.profiles.length > 1) d.profiles.forEach((p, i) => traces.push({ x, y: p, name: fmtE(d.results[i].energy), mode: 'lines', line: { width: 1.2, color: C(i + 1) } }));
    traces.push({ x, y: d.total, name: d.profiles.length > 1 ? 'Total' : fmtE(d.results[0].energy), mode: 'lines', line: { width: 2, color: C(0), dash: 'solid' }, fill: d.profiles.length > 1 ? null : 'tozeroy', fillcolor: rgba(C(0), 0.12) });
    const ymax = Math.max(...d.total);
    plot('plot-im-profile', traces, { title: `Implanted ${this.title}`, xaxis: { title: `Depth (${lu})` }, yaxis: { type: log ? 'log' : 'linear', title: 'Concentration (atoms/cm³)', range: log ? [Math.log10(ymax) - 4, Math.log10(ymax) + 0.3] : undefined } });
    const rows = this.curve.rows;
    const E = rows.map(r => r.energy);
    const [cu, cs] = lenUnit(Math.max(...rows.map(r => r.range)));
    plot('plot-im-curve', [
      { x: E, y: rows.map(r => r.range / cs), name: 'Projected range R<sub>p</sub>', mode: 'lines', line: { color: C(0), width: 2 } },
      { x: E.concat(E.slice().reverse()), y: rows.map(r => (r.range + r.long) / cs).concat(rows.map(r => Math.max(r.range - r.long, 0) / cs).reverse()), name: 'R<sub>p</sub> ± ΔR<sub>p</sub>', fill: 'toself', mode: 'lines', line: { width: 0 }, fillcolor: rgba(C(0), 0.13), hoverinfo: 'skip' },
      { x: E, y: rows.map(r => r.long / cs), name: 'Longitudinal straggling ΔR<sub>p</sub>', mode: 'lines', line: { color: C(1) } },
      { x: E, y: rows.map(r => r.lat / cs), name: 'Lateral straggling', mode: 'lines', line: { color: C(2) } },
      { x: d.results.map(r => r.energy), y: d.results.map(r => r.range / cs), name: 'Implant energies', mode: 'markers', marker: { size: 9, color: C(3), line: { color: cssVar('--panel'), width: 1.5 } } },
    ], { title: `Range and straggling of ${this.title}`, xaxis: { type: 'log', title: 'Ion energy (keV)' }, yaxis: { type: 'log', title: `Distance (${cu})` } });
    const cards = $('#im-cards'); cards.innerHTML = '';
    const peak = Math.max(...d.total);
    const peakx = d.depth[d.total.indexOf(peak)];
    cards.append(card('Peak concentration', fmt(peak, 3) + ' cm⁻³', 'at ' + fmtLen(peakx)),
      card('Peak atomic fraction', fmt(100 * peak / d.target_atom_density, 3) + ' %', 'of ' + fmt(d.target_atom_density, 3) + ' atoms/cm³'),
      card('Total dose', fmt(d.results.reduce((a, r) => a + r.dose, 0), 3) + ' cm⁻²'));
    const t = $('#im-tab'); t.innerHTML = '';
    t.append(h('tr', null, ...['Energy', 'Dose (cm⁻²)', 'Rp', 'ΔRp (long.)', 'Lateral', 'Se (MeV cm²/mg)', 'Sn (MeV cm²/mg)'].map(x => h('th', null, x))));
    for (const r of d.results) t.append(h('tr', null, h('td', null, fmtE(r.energy)), h('td', null, fmt(r.dose, 3)), h('td', null, fmtLen(r.range)), h('td', null, fmtLen(r.long)), h('td', null, fmtLen(r.lat)), h('td', null, fmt(r.se)), h('td', null, fmt(r.sn))));
  },
};
function card(k, v, s) { return h('div', { class: 'card' }, h('div', { class: 'k' }, k), h('div', { class: 'v' }, v), s ? h('div', { class: 's' }, s) : null); }

// ----------------------------------------------------------------- TRIM
const Trim = {
  layers: [],
  init() {
    this.ion = new IonPicker($('[data-ion="tr"]'), { z: 18, trim: true, energy: 50 });
    this.addLayer({ z: 14, width: 2000 });
    $('#tr-addlayer').onclick = () => this.addLayer({ z: 14, width: 1000 });
    $('#tr-run').onclick = () => this.run();
    $('#tr-stop').onclick = () => this.job && api(`/api/trim/${this.job}/cancel`, {});
    ['tr-view', 'tr-showrec', 'tr-showstop'].forEach(id => $('#' + id).addEventListener('change', () => this.draw()));
    $('#tr-png').onclick = () => { const a = h('a', { href: $('#tr-canvas').toDataURL('image/png'), download: 'trajectories.png' }); a.click(); };
    window.addEventListener('resize', () => this.draw());
    setEmpty('tr', true);
  },
  addLayer(o) {
    const wrap = h('div', { class: 'layer' });
    const L = { wrap };
    const name = h('input', { class: 'lname', value: '', placeholder: 'Layer name' });
    const width = h('input', { type: 'number', step: 'any', min: 0, value: o.width >= 1e4 ? o.width / 1e4 : o.width / 10 });
    const wu = h('select', null, ...[['Å', 1], ['nm', 10], ['µm', 1e4], ['mm', 1e7]].map(([u, f]) => h('option', { value: f, selected: (o.width >= 1e4 ? 1e4 : 10) === f }, u)));
    const tools = h('div', { class: 'tools' },
      h('button', { class: 'x-btn', title: 'Move up', onclick: () => this.move(L, -1) }, '↑'),
      h('button', { class: 'x-btn', title: 'Move down', onclick: () => this.move(L, 1) }, '↓'),
      h('button', { class: 'x-btn', title: 'Remove layer', onclick: () => { if (this.layers.length > 1) { this.layers.splice(this.layers.indexOf(L), 1); wrap.remove(); this.renumber(); } } }, '✕'));
    const badge = h('span', { class: 'badge' });
    wrap.append(h('div', { class: 'layer-head' }, badge, name, tools));
    wrap.append(h('div', { class: 'grid2' }, h('label', null, 'Width', h('span', { class: 'unit-input' }, width, wu))));
    const ed = h('div');
    wrap.append(ed);
    L.tgt = new TargetEditor(ed, { z: o.z, trim: true });
    L.tgt.onchange = () => { if (!L.userName) name.value = L.tgt.label(); };
    name.value = L.tgt.label();
    name.addEventListener('input', () => { L.userName = true; });
    L.get = () => {
      const t = L.tgt.get();
      t.name = name.value || t.name;
      t.width = (parseFloat(width.value) || 0) * parseFloat(wu.value);
      return t;
    };
    L.badge = badge;
    this.layers.push(L);
    $('#tr-layers').append(wrap);
    this.renumber();
  },
  move(L, d) {
    const i = this.layers.indexOf(L), j = i + d;
    if (j < 0 || j >= this.layers.length) return;
    this.layers.splice(i, 1); this.layers.splice(j, 0, L);
    const c = $('#tr-layers'); this.layers.forEach(l => c.append(l.wrap));
    this.renumber();
  },
  renumber() { this.layers.forEach((l, i) => { l.badge.textContent = i + 1; }); },
  async run() {
    let q;
    try {
      const ion = this.ion.get();
      ion.energy = this.ion.energy;
      q = {
        ion, layers: this.layers.map(l => l.get()),
        number_ions: Math.max(1, parseInt($('#tr-n').value) || 1), calculation: +$('#tr-calc').value,
        seed: $('#tr-seed').value === '' ? null : parseInt($('#tr-seed').value),
        trajectories: Math.max(0, Math.min(500, parseInt($('#tr-traj').value) || 0)),
        plot_xmin: (parseFloat($('#tr-xmin').value) || 0) * 10,
        plot_xmax: (parseFloat($('#tr-xmax').value) || 0) * 10,
        exyz: Math.max(0, parseFloat($('#tr-exyz').value) || 0),
        collisions: +$('#tr-coll').value,
        title: `${this.ion.label()} ${fmtE(ion.energy)} into ${this.layers.map(l => l.tgt.label()).join('/')}`,
      };
      if (!(ion.energy > 0)) throw new Error('enter the ion energy');
    } catch (e) { toast(e.message, true); return; }
    const btn = $('#tr-run');
    busy(btn, true, 'Running…');
    try {
      const j = await api('/api/trim', q);
      this.job = j.id; this.q = q; this.result = null; this.live = null;
      setEmpty('tr', false);
      $('#tr-progress').classList.remove('hidden');
      $('#tr-stop').classList.remove('hidden');
      this.poll();
    } catch (e) { toast(e.message, true); busy(btn, false); }
  },
  async poll() {
    const jid = this.job;
    let s;
    try { s = await api(`/api/trim/${jid}?live=1`); } catch (e) { toast(e.message, true); busy($('#tr-run'), false); return; }
    if (jid !== this.job) return;
    const pct = 100 * s.done / s.total;
    $('#tr-progress .fill').style.width = pct.toFixed(1) + '%';
    const rate = s.elapsed > 0 ? s.done / s.elapsed : 0;
    const eta = rate > 0 ? (s.total - s.done) / rate : 0;
    $('#tr-progress .ptext').textContent = `${s.done} / ${s.total} ions · ${s.elapsed.toFixed(1)} s` + (s.status === 'running' && rate ? ` · ${rate.toFixed(1)} ions/s · ${eta.toFixed(0)} s left` : '') + (s.live ? ` · B ${s.live.backscattered} · T ${s.live.transmitted}` : '');
    if (s.status === 'finished' || s.status === 'cancelled' || s.status === 'error') $('#tr-stop').classList.add('hidden');
    if (s.status === 'finished' || s.status === 'cancelled') {
      busy($('#tr-run'), false);
      $('#tr-progress .ptext').textContent = `${s.result.ions} ions in ${s.elapsed.toFixed(1)} s` + (s.status === 'cancelled' ? ' (stopped)' : '');
      this.result = s.result;
      this.render();
      return;
    }
    if (s.status === 'error') { busy($('#tr-run'), false); toast('TRIM failed: ' + s.error, true); return; }
    if (s.live) { this.live = s.live; this.draw(); }
    setTimeout(() => this.poll(), 700);
  },
  geometry() {
    const r = this.result;
    const layers = r ? r.layers : (() => { let t = 0; return this.q.layers.map(l => { const o = { name: l.name, top: t, bottom: t + l.width }; t += l.width; return o; }); })();
    return layers;
  },
  draw() {
    const cv = $('#tr-canvas');
    if (!cv || !cv.offsetParent || !this.q) return;
    const dpr = window.devicePixelRatio || 1;
    const W = cv.clientWidth, H = cv.clientHeight;
    cv.width = W * dpr; cv.height = H * dpr;
    const g = cv.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const src = this.result || this.live || {};
    const trajs = src.trajectories || [];
    const view = $('#tr-view').value;
    const layers = this.geometry();
    const depth = layers[layers.length - 1].bottom;
    const ia = view === 'yz' ? 1 : 0, ib = view === 'xz' ? 2 : (view === 'yz' ? 2 : 1);
    // extent
    let amin = 0, amax = depth, bmax = 1;
    const consider = (p) => { amin = Math.min(amin, p[ia]); amax = Math.max(amax, p[ia]); bmax = Math.max(bmax, Math.abs(p[ib])); };
    trajs.forEach(t => t.points.forEach(consider));
    if (view === 'yz') { amin = -bmax; amax = bmax; }
    else {
      // zoom on the region the ions reach (at least the first layer)
      let reach = layers[0].bottom;
      trajs.forEach(t => t.points.forEach(p => { reach = Math.max(reach, p[0]); }));
      if (this.result && this.result.stopped) reach = Math.max(reach, ...this.result.stopped.x.slice(0, 5000));
      amin = Math.max(amin, -0.25 * depth); amax = Math.min(Math.max(reach * 1.08, layers[0].bottom), depth * 1.15 + 1);
    }
    bmax = Math.max(bmax, (amax - amin) * 0.15);
    const narrow = isNarrow();
    const pad = narrow ? { l: 46, r: 8, t: 8, b: 38 } : { l: 64, r: 16, t: 16, b: 44 };
    const fs = narrow ? 11 : 12.5;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    let sx = pw / (amax - amin), sy = ph / (2 * bmax);
    if (view === 'yz') { const s = Math.min(sx, sy); sx = sy = s; }
    const X = (a) => pad.l + (a - amin) * sx, Y = (b) => pad.t + ph / 2 - b * sy;
    const text = cssVar('--text'), muted = cssVar('--muted'), line = cssVar('--line');
    // layers
    if (view !== 'yz') {
      // layers: alternate light grey bands, boundaries dotted
      layers.forEach((l, i) => {
        const x0 = Math.max(X(l.top), pad.l), x1 = Math.min(X(l.bottom), W - pad.r);
        if (x1 <= x0) return;
        g.fillStyle = i % 2 ? cssVar('--grid') : rgba(cssVar('--s1'), 0.05); g.fillRect(x0, pad.t, x1 - x0, ph);
        g.strokeStyle = muted; g.setLineDash([2, 3]); g.beginPath(); g.moveTo(X(l.bottom), pad.t); g.lineTo(X(l.bottom), pad.t + ph); g.stroke(); g.setLineDash([]);
        g.fillStyle = muted; g.font = `${fs}px ${FONT}`; g.textAlign = 'left'; g.fillText(l.name, x0 + 6, pad.t + 16);
      });
      g.strokeStyle = text; g.beginPath(); g.moveTo(X(0), pad.t); g.lineTo(X(0), pad.t + ph); g.stroke();
    }
    // axes
    const [lu, ls] = lenUnit(Math.max(amax - amin, 2 * bmax));
    // frame with inside ticks on all four sides, labels left and bottom
    g.fillStyle = text; g.font = `${fs}px ${FONT}`; g.strokeStyle = text; g.lineWidth = 1;
    const ticks = (lo, hi, n) => { const span = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(span))); const step = [1, 2, 5, 10].map(m => m * p).find(s => s >= span) || span; const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) out.push(Math.abs(v) < 1e-9 * step ? 0 : v); return out; };
    const tl = 6;
    g.textAlign = 'center'; g.textBaseline = 'top';
    for (const v of ticks(amin, amax, narrow ? 5 : 8)) { const x = Math.round(X(v)) + 0.5; g.fillText(fmt(v / ls, 3), x, pad.t + ph + 5); g.beginPath(); g.moveTo(x, pad.t + ph); g.lineTo(x, pad.t + ph - tl); g.moveTo(x, pad.t); g.lineTo(x, pad.t + tl); g.stroke(); }
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (const v of ticks(-bmax, bmax, narrow ? 5 : 6)) { const y = Math.round(Y(v)) + 0.5; g.fillText(fmt(v / ls, 3), pad.l - 5, y); g.beginPath(); g.moveTo(pad.l, y); g.lineTo(pad.l + tl, y); g.moveTo(pad.l + pw, y); g.lineTo(pad.l + pw - tl, y); g.stroke(); }
    g.strokeRect(pad.l + 0.5, pad.t + 0.5, pw - 1, ph - 1);
    const an = { xy: ['Depth', 'Y'], xz: ['Depth', 'Z'], yz: ['Y', 'Z'] }[view];
    g.textAlign = 'center'; g.textBaseline = 'alphabetic';
    g.fillText(`${an[0]} (${lu})`, pad.l + pw / 2, H - 6);
    g.save(); g.translate(narrow ? 11 : 14, pad.t + ph / 2); g.rotate(-Math.PI / 2); g.fillText(`${an[1]} (${lu})`, 0, 0); g.restore();
    g.textAlign = 'left';
    g.save(); g.beginPath(); g.rect(pad.l, pad.t, pw, ph); g.clip();
    // recoils
    if ($('#tr-showrec').checked && this.result && this.result.recoil_segments) {
      g.strokeStyle = rgba(cssVar('--s2'), 0.55); g.lineWidth = 0.7; g.beginPath();
      for (const s of this.result.recoil_segments) { g.moveTo(X(s[ia]), Y(s[ib])); g.lineTo(X(s[3 + ia]), Y(s[3 + ib])); }
      g.stroke();
    }
    // ion paths
    const col = { S: cssVar('--s1'), T: cssVar('--s3'), B: cssVar('--s4') };
    g.lineWidth = 1;
    for (const t of trajs) {
      g.strokeStyle = col[t.fate] || col.S; g.globalAlpha = 0.8; g.beginPath();
      t.points.forEach((p, i) => { const x = X(p[ia]), y = Y(p[ib]); if (i) g.lineTo(x, y); else g.moveTo(x, y); });
      g.stroke();
    }
    g.globalAlpha = 1;
    // final positions
    if ($('#tr-showstop').checked && this.result && this.result.stopped) {
      const st = this.result.stopped, n = Math.min(st.x.length, 8000);
      g.fillStyle = cssVar('--text'); g.globalAlpha = 0.55;
      const P = [st.x, st.y, st.z];
      for (let i = 0; i < n; i++) g.fillRect(X(P[ia][i]) - 1, Y(P[ib][i]) - 1, 2, 2);
      g.globalAlpha = 1;
    }
    g.restore();
    $('#tr-legend').innerHTML = `<div><span style="background:${col.S}"></span>stopped ions</div><div><span style="background:${col.T}"></span>transmitted</div><div><span style="background:${col.B}"></span>backscattered</div>` +
      ($('#tr-showrec').checked && this.result && this.result.recoil_segments && this.result.recoil_segments.length ? `<div><span style="background:${cssVar('--s2')}"></span>recoil atoms</div>` : '') +
      `<div class="muted">${trajs.length} trajectories${this.result ? '' : ' (live)'}</div>`;
  },
  render() {
    const r = this.result;
    this.draw();
    const sm = r.summary;
    const el = (z) => S.byZ[z].symbol;
    const elems = r.elements.map(e => el(e.z));
    const lay = r.layers;
    const depthMax = lay[lay.length - 1].bottom;
    const [lu, ls] = lenUnit(depthMax);
    const shapes = lay.slice(0, -1).map(l => ({ type: 'line', x0: l.bottom / ls, x1: l.bottom / ls, yref: 'paper', y0: 0, y1: 1, line: { color: cssVar('--muted'), dash: 'dot', width: 1 } }));
    const X = (d) => d.map(v => v / ls);
    // zoom on the depths where something happens
    let reach = lay[0].bottom;
    const scan = (dep, ...cols) => { for (let i = dep.length - 1; i >= 0; i--) if (cols.some(c => c && c[i] > 0)) { reach = Math.max(reach, dep[i]); return; } };
    if (r.range) scan(r.range.depth, r.range.ions, ...r.range.recoils);
    if (r.vacancy) scan(r.vacancy.depth, r.vacancy.knock_ons, ...r.vacancy.vacancies);
    if (r.ioniz) scan(r.ioniz.depth, r.ioniz.ions, r.ioniz.recoils);
    const binw = r.range ? (r.range.depth[1] - r.range.depth[0]) : 0;
    const xa = this.q.plot_xmax > 0 ? { title: `Depth (${lu})` } : { title: `Depth (${lu})`, range: [0, Math.min(depthMax, reach + 2 * binw) / ls] };
    // ranges
    $('#tr-range-cards').innerHTML = '';
    const LR = (v) => fmtLenAs(v, sm.range);
    $('#tr-range-cards').append(card('Mean projected range', LR(sm.range)), card('Straggling', LR(sm.straggling)),
      card('Lateral range', LR(sm.lateral), 'straggling ' + LR(sm.lateral_straggling)), card('Radial range', LR(sm.radial), 'straggling ' + LR(sm.radial_straggling)),
      card('Ions', r.ions, `${sm.backscattered} backscattered · ${sm.transmitted} transmitted`));
    const nostop = r.ions > 0 && sm.backscattered + sm.transmitted >= r.ions;
    if (nostop) ['plot-tr-range', 'plot-tr-recoil', 'plot-tr-lat', 'plot-tr-scatter'].forEach((id, i) => clearPlot(id, i ? '' : '<p class="hint">No ion stopped in the target: all were transmitted or backscattered.</p>'));
    if (r.range && !nostop) {
      const g = r.range.depth.map(x => (sm.straggling > 0 ? 1e8 * (r.ions - sm.backscattered - sm.transmitted) / r.ions / (Math.sqrt(2 * Math.PI) * sm.straggling) * Math.exp(-0.5 * ((x - sm.range) / sm.straggling) ** 2) : 0));
      plot('plot-tr-range', [
        stepTrace(X(r.range.depth), r.range.ions, 'Ions (TRIM)', C(0), { line: { shape: 'hvh', color: C(0), width: 1.5, dash: 'solid' } }),
        { x: X(r.range.depth), y: g, name: 'Gaussian (R<sub>p</sub>, ΔR<sub>p</sub>)', mode: 'lines', line: { color: C(1), width: 1.6, dash: 'dash' } },
      ], { title: 'Ion range distribution', xaxis: xa, yaxis: { title: '(atoms/cm³) / (atoms/cm²)', rangemode: 'tozero' }, shapes });
      if (r.range.recoils.length && r.range.recoils.some(c => c.some(v => v > 0))) {
        plot('plot-tr-recoil', r.range.recoils.map((c, i) => ({ x: X(r.range.depth), y: c, name: elems[i] + ' recoils', mode: 'lines', line: { color: C(i + 1) } })),
          { title: 'Final distribution of recoil atoms (target atoms displaced)', xaxis: xa, yaxis: { title: '(atoms/cm³) / (atoms/cm²)' }, shapes });
      } else clearPlot('plot-tr-recoil', '<p class="hint">Recoil atom distributions need a full damage cascade calculation.</p>');
    }
    // damage
    $('#tr-dmg-cards').innerHTML = '';
    $('#tr-dmg-cards').append(card('Vacancies / ion', fmt(sm.vacancies, 4), this.q.calculation === 1 ? 'Kinchin-Pease estimate' : null));
    if (sm.displacements != null) $('#tr-dmg-cards').append(card('Displacements / ion', fmt(sm.displacements, 4)), card('Replacements / ion', fmt(sm.replacements, 4)));
    $('#tr-dmg-cards').append(card('Sputtered atoms / ion', fmt(sm.sputtered / r.ions, 4), sm.sputtered + ' atoms'));
    if (r.vacancy) {
      const v = r.vacancy;
      const tr = [];
      if (v.vacancies.length >= 2 && this.q.calculation === 2) {
        tr.push({ x: X(v.depth), y: v.knock_ons, name: 'Knock-ons (PKA)', mode: 'lines', line: { color: C(0) } });
        v.vacancies.forEach((c, i) => tr.push({ x: X(v.depth), y: c, name: elems[i] + ' vacancies', mode: 'lines', line: { color: C(i + 1) } }));
      } else {
        tr.push({ x: X(v.depth), y: v.knock_ons, name: 'Vacancies by ions', mode: 'lines', line: { color: C(0) } });
        if (v.vacancies[0]) tr.push({ x: X(v.depth), y: v.vacancies[0], name: 'Vacancies by recoils', mode: 'lines', line: { color: C(1) } });
      }
      plot('plot-tr-vac', tr, { title: 'Vacancy production', xaxis: xa, yaxis: { title: 'Vacancies / (Å · ion)' }, shapes });
    }
    if (r.novac && r.novac.number.length) plot('plot-tr-novac', [{ x: X(r.novac.depth), y: r.novac.number, name: 'Replacement collisions', mode: 'lines', line: { color: C(0) } }], { title: 'Replacement collisions', xaxis: xa, yaxis: { title: 'Number / (Å · ion)' }, shapes });
    else clearPlot('plot-tr-novac', '');
    // energy
    if (r.exyz && r.exyz.length) {
      plot('plot-tr-exyz', r.exyz.map((c, i) => ({ x: X(c.x), y: c.e, name: 'ion ' + c.n, mode: 'lines', showlegend: false, line: { color: C(0), width: 0.8, dash: 'solid' }, opacity: 0.55 })),
        { title: `Ion energy vs depth (EXYZ.txt, ${r.exyz.length} ions)`, xaxis: xa, yaxis: { title: 'Energy (keV)' }, shapes, hovermode: 'closest' });
    } else clearPlot('plot-tr-exyz', '');
    if (r.ioniz) plot('plot-tr-ioniz', [{ x: X(r.ioniz.depth), y: r.ioniz.ions, name: 'by ions', mode: 'lines', line: { color: C(0) } }, { x: X(r.ioniz.depth), y: r.ioniz.recoils, name: 'by recoils', mode: 'lines', line: { color: C(1) } }], { title: 'Ionization', xaxis: xa, yaxis: { title: 'eV / (Å · ion)' }, shapes });
    if (r.phonons) plot('plot-tr-phon', [{ x: X(r.phonons.depth), y: r.phonons.ions, name: 'by ions', mode: 'lines', line: { color: C(0) } }, { x: X(r.phonons.depth), y: r.phonons.recoils, name: 'by recoils', mode: 'lines', line: { color: C(1) } }], { title: 'Phonons', xaxis: xa, yaxis: { title: 'eV / (Å · ion)' }, shapes });
    if (r.e2recoil) plot('plot-tr-e2r', [{ x: X(r.e2recoil.depth), y: r.e2recoil.ions, name: 'from ions', mode: 'lines', line: { color: C(0) } }].concat(r.e2recoil.absorbed.map((c, i) => ({ x: X(r.e2recoil.depth), y: c, name: r.e2recoil.absorbed.length > 1 ? 'absorbed by ' + elems[i] : 'absorbed by recoils', mode: 'lines', line: { color: C(i + 1) } }))), { title: 'Energy transferred to recoils', xaxis: xa, yaxis: { title: 'eV / (Å · ion)' }, shapes });
    // lateral
    if (r.lateral && !nostop) {
      const l = r.lateral;
      plot('plot-tr-lat', [
        { x: X(l.depth), y: l.lateral_range.map(v => v / ls), name: 'Lateral projected range', mode: 'lines', line: { color: C(0) } },
        { x: X(l.depth), y: l.lateral_straggling.map(v => v / ls), name: 'Projected straggling', mode: 'lines', line: { color: C(1) } },
        { x: X(l.depth), y: l.radial_range.map(v => v / ls), name: 'Radial range', mode: 'lines', line: { color: C(2) } },
        { x: X(l.depth), y: l.radial_straggling.map(v => v / ls), name: 'Radial straggling', mode: 'lines', line: { color: C(3) } },
      ], { title: 'Lateral spread of the stopped ions vs depth', xaxis: xa, yaxis: { title: `(${lu})` }, shapes });
    }
    if (r.stopped && !nostop) {
      const s = r.stopped;
      // ion density seen along the beam: empty bins in the panel colour, then one hue light -> dark
      const m = Math.max(1e-9, ...s.y.map(Math.abs), ...s.z.map(Math.abs)) / ls * 1.05;
      const c0 = C(0), pc = cssVar('--panel');
      plot('plot-tr-scatter', [{ x: s.y.map(v => v / ls), y: s.z.map(v => v / ls), type: 'histogram2d', name: 'ions',
        xbins: { start: -m, end: m, size: m / 30 }, ybins: { start: -m, end: m, size: m / 30 }, zmin: 0,
        colorscale: [[0, pc], [1e-6, rgba(c0, 0.15)], [0.5, rgba(c0, 0.6)], [1, c0]], hovertemplate: 'Y %{x}<br>Z %{y}<br>%{z} ions<extra></extra>',
        colorbar: { title: { text: 'ions', side: 'right' }, thickness: 12, outlinewidth: 1, outlinecolor: cssVar('--text'), ticks: 'inside', tickcolor: cssVar('--text'), tickfont: { color: cssVar('--text') } } }],
        { title: `Final ion positions seen along the beam (${s.count} ions)`, xaxis: { title: `Y (${lu})`, range: [-m, m], constrain: 'domain', showgrid: false }, yaxis: { title: `Z (${lu})`, range: [-m, m], scaleanchor: 'x', constrain: 'domain', showgrid: false }, hovermode: 'closest' });
    }
    // transmitted / backscattered / sputtered
    const oc = $('#tr-out-cards'); oc.innerHTML = '';
    const tr = [], ang = [], outs = [];
    const add = (o, name, c) => {
      if (!o) return;
      const mean = o.energy.reduce((a, b) => a + b, 0) / o.energy.length;
      oc.append(card(name, o.count, 'mean energy ' + fmtE(mean / 1000)));
      const th = o.cx.map((c0, i) => (name === 'Backscattered' ? 180 / Math.PI * Math.acos(Math.min(1, -c0)) : 180 / Math.PI * Math.acos(Math.min(1, Math.abs(c0)))));
      outs.push({ e: o.energy.map(v => v / 1000), th, name, c });
    };
    add(r.transmitted, 'Transmitted', C(2)); add(r.backscattered, 'Backscattered', C(3)); add(r.sputtered, 'Sputtered', C(1));
    if (outs.length) {
      const emax = Math.max(...outs.map(o => Math.max(...o.e)));
      outs.forEach((o, i) => {
        const he = histogram(o.e, 0, emax * 1.0001, 60), ha = histogram(o.th, 0, 90, 45);
        const dash = DASHES[i];
        tr.push(stepTrace(he.x, he.y, o.name, o.c, { line: { shape: 'hvh', color: o.c, width: 1.5, dash } }));
        ang.push(stepTrace(ha.x, ha.y, o.name, o.c, { line: { shape: 'hvh', color: o.c, width: 1.5, dash } }));
      });
    }
    if (tr.length) {
      plot('plot-tr-espec', tr, { title: 'Energy spectra of atoms leaving the target', xaxis: { title: 'Energy (keV)', rangemode: 'tozero' }, yaxis: { title: 'Counts', rangemode: 'tozero' } });
      plot('plot-tr-ang', ang, { title: 'Angle to the surface normal (exit direction)', xaxis: { title: 'Angle (deg)', range: [0, 90], dtick: 15 }, yaxis: { title: 'Counts', rangemode: 'tozero' } });
    } else {
      oc.append(card('No particles left the target', '0'));
      clearPlot('plot-tr-espec', ''); clearPlot('plot-tr-ang', '');
    }
    // summary & files
    const sc = $('#tr-sum-cards'); sc.innerHTML = '';
    sc.append(card('Ions', r.ions), card('Projected range', fmtLen(sm.range), 'straggling ' + fmtLenAs(sm.straggling, sm.range)), card('Vacancies / ion', fmt(sm.vacancies, 4)), card('Backscattered', sm.backscattered), card('Transmitted', sm.transmitted), card('Sputtered', sm.sputtered));
    const job = this.job;
    if (LS) {
      $('#tr-zip').href = '#';
      $('#tr-zip').onclick = async (ev) => {
        ev.preventDefault();
        const r = await LS.send({ kind: 'zip', job });
        if (r.status !== 200) { toast(JSON.parse(r.text).detail, true); return; }
        const a = h('a', { href: URL.createObjectURL(new Blob([r.bytes], { type: 'application/zip' })), download: `lunasrim_${job}.zip` });
        a.click();
      };
    } else $('#tr-zip').href = `/api/trim/${job}/files.zip`;
    const ul = $('#tr-files'); ul.innerHTML = '';
    const fileText = async (f) => LS ? (await LS.send({ kind: 'file', job, name: f })).text : (await fetch(`/api/trim/${job}/file/${f}`)).text();
    for (const f of r.files) ul.append(h('li', null, h('a', { href: '#', onclick: async (ev) => { ev.preventDefault(); $('#tr-filepre').textContent = await fileText(f); } }, f)));
  },
};

// -------------------------------------------------------------- compare
const CMP = {
  curves: [],
  init() {
    this.ion = new IonPicker($('[data-ion="cmp"]'), { z: 2 });
    this.tgt = new TargetEditor($('[data-target="cmp"]'), { z: 79 });
    $('#cmp-add').onclick = () => this.add();
    $('#cmp-q').onchange = () => this.render();
    $('#cmp-peru').onchange = () => this.render();
    setEmpty('cmp', true);
  },
  async add() {
    const btn = $('#cmp-add');
    let q;
    try {
      q = { ion: this.ion.get(), target: this.tgt.get(), emin: parseFloat($('#cmp-emin').value), emax: parseFloat($('#cmp-emax').value), units: 5 };
    } catch (e) { toast(e.message, true); return; }
    busy(btn, true);
    try {
      const d = await api('/api/stopping', q);
      this.curves.push({ label: `${this.ion.label()} → ${this.tgt.label()}`, mass: q.ion.mass, rows: d.rows, slot: this.freeSlot() });
      setEmpty('cmp', false);
      this.render();
    } catch (e) { toast(e.message, true); } finally { busy(btn, false); }
  },
  freeSlot() {
    for (let i = 0; ; i++) if (!this.curves.some(c => c.slot === i)) return i;
  },
  render() {
    const ul = $('#cmp-list'); ul.innerHTML = '';
    this.curves.forEach((c, i) => ul.append(h('li', null, h('span', { class: 'sw', style: `background:${C(c.slot)}` }), h('span', { class: 'lbl' }, c.label), h('button', { class: 'x-btn', onclick: () => { this.curves.splice(i, 1); this.render(); } }, '✕'))));
    if (!this.curves.length) { setEmpty('cmp', true); return; }
    const q = $('#cmp-q').value, peru = $('#cmp-peru').checked;
    const val = { total: r => r.se + r.sn, se: r => r.se, sn: r => r.sn, range: r => r.range * 1e-4, long: r => r.long * 1e-4, lat: r => r.lat * 1e-4, rel: r => r.long / r.range }[q];
    const ytitle = { total: 'Total stopping (MeV cm²/mg)', se: 'Electronic stopping (MeV cm²/mg)', sn: 'Nuclear stopping (MeV cm²/mg)', range: 'Projected range (µm)', long: 'Longitudinal straggling (µm)', lat: 'Lateral straggling (µm)', rel: 'ΔRp / Rp' }[q];
    plot('plot-cmp', this.curves.map(c => ({ x: c.rows.map(r => peru ? r.energy / c.mass : r.energy), y: c.rows.map(val), name: c.label, mode: 'lines', line: { color: C(c.slot), dash: DASHES[c.slot % DASHES.length] } })),
      { xaxis: { type: 'log', title: peru ? 'Energy (keV/u)' : 'Energy (keV)' }, yaxis: { type: q === 'rel' ? 'linear' : 'log', title: ytitle } });
  },
};

// ------------------------------------------------------------------ theme
function initTheme() {
  const saved = (() => { try { return localStorage.getItem('lunasrim-theme'); } catch (e) { return null; } })();
  if (saved) document.documentElement.dataset.theme = saved;
  $('#themeToggle').onclick = () => {
    const cur = document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const nxt = cur === 'dark' ? 'light' : 'dark';
    document.documentElement.dataset.theme = nxt;
    try { localStorage.setItem('lunasrim-theme', nxt); } catch (e) { /* private mode */ }
    if (SR.data) SR.render(); if (IM.data) IM.render(); if (Trim.result) Trim.render(); else Trim.draw(); CMP.render();
  };
}

// ------------------------------------------------------------------- boot
async function boot() {
  initTheme();
  const loading = $('#loading');
  if (LS && loading) loading.classList.remove('hidden');
  try {
    const [els, comps, info] = await Promise.all([api('/api/elements'), api('/api/compounds'), api('/api/info')]);
    S.elements = els; S.compounds = comps; S.info = info;
    for (const e of els) { S.byZ[e.z] = e; S.bySym[e.symbol] = e; }
  } catch (e) {
    toast((LS ? 'Cannot start the Python engine: ' : 'Cannot reach the LunaSRIM server: ') + e.message, true);
    if (loading) loading.textContent = 'The Python engine could not be started: ' + e.message;
    return;
  }
  if (loading) loading.classList.add('hidden');
  Picker.build();
  SR.init(); IM.init(); Trim.init(); CMP.init();
  initTabs();
  $('#about-version').textContent = S.info.static
    ? `lunasrim ${S.info.version} · running in your browser (Pyodide), no data leaves your computer`
    : `lunasrim ${S.info.version} · up to ${S.info.max_ions} ions per TRIM calculation on this server`;
  $('#tr-n').max = S.info.max_ions;
}
document.addEventListener('DOMContentLoaded', boot);
