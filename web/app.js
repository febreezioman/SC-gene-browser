// Cell Read Browser: build a cell selection on the UMAP, then load those cells'
// reads in igv.js.
//
// One plotly panel per condition (clusters are numbered per condition, so they
// are never shown on a shared colour key). Selecting is cheap and local; reads
// load only when "Load reads" is pressed. igv.js then gets one alignment track
// per condition whose filter keeps reads whose CB tag is a selected cell. The
// filter is a function, which igv.js applies before downsampling.

import igv from "https://cdn.jsdelivr.net/npm/igv@3.8.8/dist/igv.esm.min.js";

const params = new URLSearchParams(location.search);
const LINE = params.get("line") || "AP2XI-2_AP2XII-1";
const DATA = `../data/${LINE}/`;
const REF = "../data/ref/";
// AP2XI-2 (TGME49_310900): a knockdown target, and light enough to load fast.
// Very highly expressed genes (SRS29B/SAG1: ~290k reads deep) take a while,
// because every read in the window is downloaded and parsed before filtering.
const START_LOCUS = params.get("locus") || "TGME49_chrXI:1,645,000-1,654,500";

// Categorical slots in fixed order (dataviz reference palette, light mode).
const SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"];
const OTHER = "#8a8984";
const DIM = "rgba(190, 189, 184, 0.35)";
// Sequential blue ramp, light -> dark.
const RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];

// Fixed category orders so a colour always means the same entity.
const ORDERS = {
  phase: ["G1a", "G1b", "S", "M", "C"],
  brady_ees: ["Tachyzoite", "Bradyzoite", "Bradyzoite/EES1-3", "EES1-3"],
  stage10: ["Tachyzoites", "BRP1+", "Merozoite A", "Merozoite B", "Merozoite C",
            "Trans. Merozoites", "Pre-gamete", "Female", "Male-1", "Male-2"],
};
const CONTINUOUS = { pseudotime: "Pseudotime", umis: "UMIs per cell" };

const state = {
  cells: null,          // columnar table from cells.json
  n: 0,
  selected: new Set(),  // global cell indices, the selection being built
  loaded: null,         // Set of cell indices whose reads are shown; null = none loaded
  loadedText: "",
  pick: "group",
  mode: "add",
  colorBy: "cluster",
  groupBy: "none",
  readColor: "strand",
  displayMode: "EXPANDED",
  trackHeight: null,    // px per read track; null = fit to the panel
  highlight: null,      // cell index of the last clicked read
};
let browser;
let readTracks = [];
let cbIndex;            // per condition: Map CB -> cell index
let hoverText;          // per cell tooltip
const plots = [];       // per condition: {div, idx: [cell indices]}

const $ = id => document.getElementById(id);
const fmt = n => n.toLocaleString();

// ---------- data and groups ----------

async function loadJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json();
}

// The variable that defines a "group": the colour variable, or the cluster
// when colouring by a continuous value.
const groupVar = () => (state.colorBy in CONTINUOUS ? "cluster" : state.colorBy);

function groupOf(i, v = groupVar()) {
  const c = state.cells;
  return v === "cluster" ? String(c.cluster[i]) : c[v][i];
}

function groupName(key, v = groupVar()) {
  return v === "cluster" ? `Cluster ${key}` : (key || "unassigned");
}

function groupOrder(cond, v = groupVar()) {
  if (v === "cluster") {
    const ids = new Set(plots[cond].idx.map(i => state.cells.cluster[i]));
    return [...ids].sort((a, b) => a - b).map(String);
  }
  return ORDERS[v];
}

function colorFor(i, cond) {
  if (state.colorBy in CONTINUOUS) {
    const [lo, hi] = continuousRange();
    const v = state.cells[state.colorBy][i];
    if (v < 0) return OTHER;
    const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo || 1)));
    return RAMP[Math.round(t * (RAMP.length - 1))];
  }
  const k = groupOrder(cond).indexOf(groupOf(i));
  return k >= 0 && k < SLOTS.length ? SLOTS[k] : OTHER;
}

const rangeCache = {};
function continuousRange() {
  const key = state.colorBy;
  if (!rangeCache[key]) {
    const v = state.cells[key].filter(x => x >= 0).sort((a, b) => a - b);
    // 1st-99th percentile so a few extreme cells do not flatten the ramp
    rangeCache[key] = [v[Math.floor(v.length * 0.01)], v[Math.floor(v.length * 0.99)]];
  }
  return rangeCache[key];
}

function groupMembers(cond, key) {
  return plots[cond].idx.filter(i => groupOf(i) === key);
}

// ---------- UMAP panels ----------

function buildPanels() {
  const c = state.cells;
  const pad = 0.5;
  const xr = [Math.min(...c.x) - pad, Math.max(...c.x) + pad];
  const yr = [Math.min(...c.y) - pad, Math.max(...c.y) + pad];
  const host = $("panels");

  c.conditions.forEach((cond, k) => {
    const idx = [];
    for (let i = 0; i < state.n; i++) if (c.cond[i] === k) idx.push(i);
    const wrap = document.createElement("div");
    wrap.className = "panel";
    wrap.innerHTML = `<h2>${cond} <span class="n">· ${fmt(idx.length)} cells</span></h2><div class="plot"></div>`;
    host.appendChild(wrap);
    const div = wrap.querySelector(".plot");
    plots.push({ div, idx });

    const axis = { visible: false };
    const flat = { selected: { marker: { opacity: 1 } }, unselected: { marker: { opacity: 1 } } };
    Plotly.newPlot(div, [
      { type: "scattergl", mode: "markers", x: [], y: [], customdata: [], text: [],
        hovertemplate: "%{text}<extra></extra>", marker: { size: 3, color: [] }, ...flat },
      { type: "scattergl", mode: "markers", x: [], y: [], customdata: [], text: [],
        hovertemplate: "%{text}<extra></extra>", marker: { size: 4.5, color: [] }, ...flat },
      // Ring for the cell of the last clicked read. WebGL, because plotly
      // draws WebGL traces above SVG ones.
      { type: "scattergl", mode: "markers", x: [], y: [], hoverinfo: "skip",
        marker: { size: 16, symbol: "circle-open", color: "#0b0b0b", line: { width: 2.5 } } },
    ], {
      margin: { l: 0, r: 0, t: 0, b: 0 }, showlegend: false, dragmode: "lasso",
      xaxis: { ...axis, range: xr }, yaxis: { ...axis, range: yr, scaleanchor: "x" },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      hovermode: "closest", hoverlabel: { bgcolor: "#fff", bordercolor: "#dcdbd6", font: { color: "#0b0b0b" } },
    }, {
      displaylogo: false, responsive: true,
      modeBarButtonsToRemove: ["toImage", "autoScale2d", "zoomIn2d", "zoomOut2d", "select2d"],
    });

    div.on("plotly_selected", ev => {
      if (!ev || !ev.points || !ev.points.length) return;
      const ids = ev.points.map(p => p.customdata).filter(i => i !== undefined);
      if (state.mode === "add") ids.forEach(i => state.selected.add(i));
      else ids.forEach(i => state.selected.delete(i));
      Plotly.relayout(div, { selections: [] });
      onSelectionChanged();
    });
    div.on("plotly_click", ev => {
      const p = ev.points && ev.points[0];
      if (!p || p.customdata === undefined) return;
      if (state.pick === "group") toggleGroup(k, groupOf(p.customdata));
      else toggleCells([p.customdata]);
    });
  });
}

function redrawPanels() {
  const c = state.cells;
  const has = state.selected.size > 0;
  plots.forEach(({ div, idx }, cond) => {
    const bg = { x: [], y: [], cd: [], col: [], tx: [] };
    const fg = { x: [], y: [], cd: [], col: [], tx: [] };
    for (const i of idx) {
      const on = !has || state.selected.has(i);
      const t = on && has ? fg : bg;
      t.x.push(c.x[i]); t.y.push(c.y[i]); t.cd.push(i); t.tx.push(hoverText[i]);
      t.col.push(on ? colorFor(i, cond) : DIM);
    }
    const h = state.highlight !== null && c.cond[state.highlight] === cond ? state.highlight : null;
    Plotly.restyle(div, {
      x: [bg.x, fg.x, h === null ? [] : [c.x[h]]],
      y: [bg.y, fg.y, h === null ? [] : [c.y[h]]],
      customdata: [bg.cd, fg.cd, []],
      text: [bg.tx, fg.tx, []],
      "marker.color": [bg.col, fg.col, "#0b0b0b"],
      selectedpoints: [null, null, null],
    });
  });
}

// ---------- legend (doubles as group selector) ----------

function renderLegend() {
  const host = $("legend");
  host.innerHTML = "";
  const c = state.cells;
  plots.forEach(({ idx }, cond) => {
    const col = document.createElement("div");
    col.className = "col";
    col.innerHTML = `<h3>${c.conditions[cond]}</h3>`;
    if (state.colorBy in CONTINUOUS) {
      const [lo, hi] = continuousRange();
      const f = v => (state.colorBy === "umis" ? fmt(Math.round(v)) : v.toFixed(2));
      col.innerHTML += `<div class="ramp" style="background:linear-gradient(90deg,${RAMP.join(",")})"></div>
        <div class="ramp-labels"><span>${f(lo)}</span><span>${CONTINUOUS[state.colorBy]}</span><span>${f(hi)}</span></div>
        <h3>Groups (clusters)</h3>`;
    }
    const counts = new Map(), sel = new Map();
    for (const i of idx) {
      const g = groupOf(i);
      counts.set(g, (counts.get(g) || 0) + 1);
      if (state.selected.has(i)) sel.set(g, (sel.get(g) || 0) + 1);
    }
    const order = groupOrder(cond).filter(k => counts.has(k));
    for (const k of counts.keys()) if (!order.includes(k)) order.push(k);
    for (const k of order) {
      const slot = groupOrder(cond).indexOf(k);
      const sw = state.colorBy in CONTINUOUS ? "transparent"
        : slot >= 0 && slot < SLOTS.length ? SLOTS[slot] : OTHER;
      const n = counts.get(k), s = sel.get(k) || 0;
      const b = document.createElement("button");
      b.className = "chip" + (s === n ? " full" : s > 0 ? " partial" : "");
      b.innerHTML = `<span class="sw" style="background:${sw}"></span>${groupName(k)}` +
        `<span class="n">${s && s < n ? `<span class="sel">${fmt(s)}</span> / ` : ""}${fmt(n)}</span>`;
      b.title = `${s === n ? "Remove" : "Add"} ${groupName(k)} (${c.conditions[cond]})`;
      b.setAttribute("aria-pressed", s === n);
      b.onclick = () => toggleGroup(cond, k);
      col.appendChild(b);
    }
    host.appendChild(col);
  });
}

// ---------- selection ----------

function toggleCells(ids) {
  const allIn = ids.every(i => state.selected.has(i));
  ids.forEach(i => (allIn ? state.selected.delete(i) : state.selected.add(i)));
  onSelectionChanged();
}

const toggleGroup = (cond, key) => toggleCells(groupMembers(cond, key));

function selectedByCondition(set) {
  const c = state.cells;
  const per = c.conditions.map(() => new Set());
  for (const i of set) per[c.cond[i]].add(c.cb[i]);
  return per;
}

// "UT: Cluster 2, Cluster 3 · IAA 24h: Cluster 4 + 12 cells"
function describeSelection(set) {
  if (set.size === 0) return { head: "All cells", groups: "" };
  const c = state.cells;
  const parts = [];
  plots.forEach(({ idx }, cond) => {
    const full = [];
    let loose = 0;
    const byGroup = new Map();
    for (const i of idx) {
      const g = groupOf(i);
      if (!byGroup.has(g)) byGroup.set(g, [0, 0]);
      const e = byGroup.get(g);
      e[0]++;
      if (set.has(i)) e[1]++;
    }
    for (const [g, [n, s]] of byGroup) {
      if (s === n) full.push(groupName(g));
      else loose += s;
    }
    if (!full.length && !loose) return;
    const bits = full.slice();
    if (loose) bits.push(`${fmt(loose)} ${full.length ? "more " : ""}cell${loose === 1 ? "" : "s"}`);
    parts.push(`${c.conditions[cond]}: ${bits.join(", ")}`);
  });
  const per = selectedByCondition(set);
  const head = `${fmt(set.size)} cells (` +
    c.conditions.map((cond, k) => `${cond} ${fmt(per[k].size)}`).join(" · ") + ")";
  return { head, groups: parts.join(" · ") };
}

function sameSet(a, b) {
  if (!a || a.size !== b.size) return false;
  for (const x of a) if (!b.has(x)) return false;
  return true;
}

function onSelectionChanged() {
  const d = describeSelection(state.selected);
  $("sel-summary").innerHTML = state.selected.size
    ? `${d.head}<span class="groups">${d.groups}</span>`
    : "No cells selected";
  redrawPanels();
  renderLegend();
  updateLoadState();
}

function updateLoadState() {
  const btn = $("load-reads");
  const current = sameSet(state.loaded, state.selected);
  btn.disabled = current;
  btn.textContent = current ? "Reads loaded"
    : state.selected.size ? `Load reads (${fmt(state.selected.size)} cells)` : "Load reads for all cells";
  const bar = $("loaded-bar");
  if (state.loaded === null) { bar.hidden = true; return; }
  bar.hidden = false;
  bar.classList.toggle("stale", !current);
  bar.innerHTML = `Showing reads for <b>${state.loadedText}</b>` +
    (current ? "" : " · selection changed: press <b>Load reads</b> to update");
}

function exportBarcodes() {
  const c = state.cells;
  const rows = ["condition\tbarcode\tcluster"];
  const ids = [...state.selected].sort((a, b) => a - b);
  ids.forEach(i => rows.push(`${c.conditions[c.cond[i]]}\t${c.cb[i]}\t${c.cluster[i]}`));
  const blob = new Blob([rows.join("\n") + "\n"], { type: "text/tab-separated-values" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `${LINE}_selection_${ids.length}cells.tsv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- igv ----------

function readTrackConfig(cond, k, cbs, all) {
  const c = state.cells;
  const label = all ? "all cells" : `${fmt(cbs.size)} cells`;
  const cfg = {
    type: "alignment", format: "bam",
    name: `${c.line} ${cond}: ${label}`,
    url: DATA + c.bams[cond], indexURL: DATA + c.bams[cond] + ".bai",
    // A function filter replaces igv's default one, so keep its default checks.
    filter: a => a.isMapped() && !a.isFailsVendorQualityCheck() && !a.isDuplicate() &&
                 (all || cbs.has(a.getTag("CB"))),
    samplingDepth: 10000,
    // igv.js 3.8.8 bug: when its block cache partly overlaps a request, it can
    // compute an empty tail (size 0) and send an open-ended "bytes=N-" Range,
    // downloading the rest of the BAM. Tracks are rebuilt per load, so the
    // cache buys nothing here anyway.
    cacheBlocks: false,
    displayMode: state.displayMode,
    showSoftClips: false, showInsertions: false,
    height: state.trackHeight || 400, maxHeight: 20000, coverageTrackHeight: 50,
    order: 10 + k,
  };
  if (state.groupBy !== "none") cfg.groupBy = state.groupBy;
  if (state.readColor !== "none") cfg.colorBy = state.readColor;
  return cfg;
}

let generation = 0;
async function loadReads() {
  if (!browser) return;
  const set = new Set(state.selected);
  state.loaded = set;
  const d = describeSelection(set);
  state.loadedText = d.groups || d.head;
  updateLoadState();
  await rebuildReadTracks();
}

async function rebuildReadTracks() {
  if (!browser || state.loaded === null) return;
  const gen = ++generation;
  for (const t of readTracks) browser.removeTrack(t);
  readTracks = [];
  const set = state.loaded;
  const per = selectedByCondition(set);
  const all = set.size === 0;
  $("reads-empty").hidden = true;
  const jobs = state.cells.conditions.map(async (cond, k) => {
    if (!all && per[k].size === 0) return; // nothing selected in this condition
    try {
      const track = await browser.loadTrack(readTrackConfig(cond, k, per[k], all));
      // A newer load started while this track was loading: drop it.
      if (gen !== generation) browser.removeTrack(track);
      else readTracks.push(track);
    } catch (err) {
      console.warn(`Read track for ${cond} failed to load`, err);
    }
  });
  await Promise.all(jobs);
  if (gen === generation) applyTrackHeight();
}

// Read tracks share the panel height unless the slider set a fixed height.
function applyTrackHeight() {
  if (!readTracks.length) return;
  let h = state.trackHeight;
  if (h === null) {
    const wrap = $("igv-wrap");
    const used = readTracks.reduce((s, t) => s + t.trackView.track.height, 0);
    const other = $("igv").scrollHeight - used;
    h = Math.max(150, Math.floor((wrap.clientHeight - other - 12) / readTracks.length));
    $("track-height").value = h;
  }
  for (const t of readTracks) t.trackView.setTrackHeight(h, true);
}

async function initIGV() {
  browser = await igv.createBrowser($("igv"), {
    reference: {
      id: "tgme49_toxodb13", name: "T. gondii ME49 (ToxoDB-13)",
      fastaURL: REF + "genome.fasta", indexURL: REF + "genome.fasta.fai",
    },
    locus: START_LOCUS,
    showNavigation: true, showSampleNames: false,
    tracks: [{
      name: "Genes (ToxoDB-13)", type: "annotation", format: "gtf",
      url: REF + "genes.gtf.gz", indexURL: REF + "genes.gtf.gz.tbi",
      displayMode: "EXPANDED", order: 1, height: 70,
    }],
  });

  // Clicking a read: find its cell and mark it on the UMAP.
  browser.on("trackclick", (track, data) => {
    if (track.type !== "alignment" || !Array.isArray(data)) return true;
    const cb = data.find(d => d && d.name === "CB");
    if (!cb) return true;
    const c = state.cells;
    const cond = c.conditions.find(x => track.config.url.endsWith(c.bams[x]));
    const i = cond === undefined ? undefined : cbIndex[c.conditions.indexOf(cond)].get(String(cb.value));
    if (i !== undefined) showReadCell(i);
    return true;
  });
}

function showReadCell(i) {
  const c = state.cells;
  state.highlight = i;
  const host = $("read-cell");
  host.hidden = false;
  host.innerHTML = `Read from cell <b>${c.cb[i]}</b> · ${c.conditions[c.cond[i]]} cluster ${c.cluster[i]}
    · ${fmt(c.umis[i])} UMIs (ringed on the UMAP) <button>Add to selection</button>`;
  host.querySelector("button").onclick = () => { state.selected.add(i); onSelectionChanged(); };
  redrawPanels();
}

// ---------- gene search ----------

let genes = [];
function initSearch() {
  const input = $("search-input");
  const list = $("search-results");
  const norm = s => s.toLowerCase().replace(/[-_\s.]/g, "");
  const index = genes.map(g => ({ g, id: norm(g[0]), syms: g[6].map(norm), desc: g[5].toLowerCase() }));
  const go = locus => { list.hidden = true; browser.search(locus); };
  const locusOf = g => {
    const pad = Math.round((g[3] - g[2]) * 0.15) + 200;
    return `${g[1]}:${Math.max(1, g[2] - pad)}-${g[3] + pad}`;
  };
  const matches = raw => {
    const q = norm(raw), ql = raw.trim().toLowerCase();
    if (q.length < 2) return [];
    const scored = [];
    for (const e of index) {
      let s = 0;
      if (e.syms.includes(q) || e.id === q) s = 4;
      else if (e.syms.some(x => x.startsWith(q))) s = 3;
      else if (e.id.includes(q)) s = 2;
      else if (ql.length >= 3 && e.desc.includes(ql)) s = 1;
      if (s) scored.push([s, e.g]);
    }
    scored.sort((a, b) => b[0] - a[0] || a[1][0].localeCompare(b[1][0]));
    return scored.slice(0, 50).map(x => x[1]);
  };
  let current = [], active = -1;
  const render = () => {
    list.innerHTML = current.map((g, j) =>
      `<li data-j="${j}" class="${j === active ? "active" : ""}">` +
      (g[6].length ? `<span class="sym">${g[6].join(" / ")}</span>` : "") +
      `<span class="gid">${g[0]}</span><span class="desc">${g[5]}</span></li>`).join("");
    list.hidden = current.length === 0;
    list.querySelectorAll("li").forEach(li => li.onclick = () => pick(current[+li.dataset.j]));
  };
  const pick = g => { input.value = g[6][0] || g[0]; go(locusOf(g)); };
  input.addEventListener("input", () => { current = matches(input.value); active = -1; render(); });
  input.addEventListener("keydown", e => {
    if (list.hidden || !current.length) return;
    if (e.key === "ArrowDown") { active = Math.min(active + 1, current.length - 1); render(); e.preventDefault(); }
    if (e.key === "ArrowUp") { active = Math.max(active - 1, 0); render(); e.preventDefault(); }
    if (e.key === "Escape") list.hidden = true;
  });
  $("search").addEventListener("submit", e => {
    e.preventDefault();
    const m = active >= 0 ? [current[active]] : matches(input.value);
    if (m.length) pick(m[0]); else go(input.value.trim());
  });
  document.addEventListener("click", e => { if (!e.target.closest("#search")) list.hidden = true; });
}

// ---------- layout: draggable splitter ----------

function initSplitter() {
  const split = $("splitter");
  const root = document.documentElement;
  try {
    const w = localStorage.getItem("cellsWidth");
    if (w) root.style.setProperty("--cells-width", w);
  } catch { /* storage unavailable */ }
  const relayout = () => {
    plots.forEach(p => Plotly.Plots.resize(p.div));
    // igv.js only re-measures on window resize; it has no public resize().
    if (browser && browser.boundWindowResizeHandler) browser.boundWindowResizeHandler();
    if (state.trackHeight === null) applyTrackHeight();
  };
  split.addEventListener("pointerdown", e => {
    split.setPointerCapture(e.pointerId);
    split.classList.add("dragging");
    const left = $("cells").getBoundingClientRect().left;
    const move = ev => {
      const w = Math.max(320, Math.min(ev.clientX - left, innerWidth - 420));
      root.style.setProperty("--cells-width", `${w}px`);
    };
    const up = () => {
      split.classList.remove("dragging");
      split.removeEventListener("pointermove", move);
      split.removeEventListener("pointerup", up);
      try { localStorage.setItem("cellsWidth", root.style.getPropertyValue("--cells-width")); } catch { }
      relayout();
    };
    split.addEventListener("pointermove", move);
    split.addEventListener("pointerup", up);
  });
  let t;
  addEventListener("resize", () => { clearTimeout(t); t = setTimeout(relayout, 150); });
}

// ---------- wiring ----------

function segmented(id, onPick) {
  const buttons = document.querySelectorAll(`#${id} button`);
  buttons.forEach(b => b.onclick = () => {
    buttons.forEach(x => x.classList.toggle("on", x === b));
    onPick(b);
  });
}

async function main() {
  const [cells, geneList] = await Promise.all([loadJSON(DATA + "cells.json"), loadJSON(REF + "genes.json")]);
  state.cells = cells;
  state.n = cells.cb.length;
  genes = geneList;
  cbIndex = cells.conditions.map(() => new Map());
  for (let i = 0; i < state.n; i++) cbIndex[cells.cond[i]].set(cells.cb[i], i);
  hoverText = cells.cb.map((_, i) =>
    `${cells.conditions[cells.cond[i]]} · Cluster ${cells.cluster[i]}<br>` +
    `${cells.phase[i]} · ${cells.stage10[i]} · ${fmt(cells.umis[i])} UMIs`);
  $("line-label").textContent = cells.line;

  buildPanels();
  onSelectionChanged();
  initSplitter();
  await initIGV();
  initSearch();

  $("color-by").onchange = e => { state.colorBy = e.target.value; onSelectionChanged(); };
  segmented("pick", b => { state.pick = b.dataset.pick; });
  segmented("mode", b => { state.mode = b.dataset.mode; });
  $("sel-clear").onclick = () => { state.selected = new Set(); onSelectionChanged(); };
  $("sel-export").onclick = exportBarcodes;
  $("load-reads").onclick = loadReads;

  $("group-by").onchange = e => { state.groupBy = e.target.value; rebuildReadTracks(); };
  $("read-color").onchange = e => { state.readColor = e.target.value; rebuildReadTracks(); };
  $("display-mode").onchange = e => {
    state.displayMode = e.target.value;
    readTracks.forEach(t => t.alignmentTrack.setDisplayMode(state.displayMode));
  };
  $("track-height").value = 400;
  $("track-height").oninput = e => { state.trackHeight = +e.target.value; applyTrackHeight(); };
  $("fit-height").onclick = () => { state.trackHeight = null; applyTrackHeight(); };
}

main().catch(err => {
  document.body.insertAdjacentHTML("afterbegin",
    `<p style="padding:12px;background:#fde8e8;margin:0">Failed to load: ${err.message}.
     Serve the project with <code>python3 scripts/serve.py</code> and run the prep scripts first.</p>`);
  throw err;
});
