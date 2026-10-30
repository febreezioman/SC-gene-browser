// Cell Read Browser: UMAP selection drives per-cell read filtering in igv.js.
//
// One plotly panel per condition (clusters are numbered per condition, so they
// are never shown on a shared colour key). igv.js gets one alignment track per
// condition whose filter keeps only reads whose CB tag is a selected cell.
// The filter is a function, which igv.js applies before downsampling.

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
  selected: new Set(),  // global cell indices
  mode: "replace",
  colorBy: "cluster",
  groupBy: "tag:CL",
  readColor: "strand",
  shift: false,
  highlight: null,      // cell index of the last clicked read
};
let browser;
let readTracks = [];
let cbIndex;            // per condition: Map CB -> cell index
const plots = [];       // per condition: {div, idx: [cell indices]}

// ---------- data ----------

async function loadJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json();
}

function categoryOf(i) {
  const c = state.cells;
  return state.colorBy === "cluster" ? String(c.cluster[i]) : c[state.colorBy][i];
}

function categoryOrder(cond) {
  if (state.colorBy === "cluster") {
    const ids = new Set(plots[cond].idx.map(i => state.cells.cluster[i]));
    return [...ids].sort((a, b) => a - b).map(String);
  }
  return ORDERS[state.colorBy];
}

function colorFor(i, cond) {
  if (state.colorBy in CONTINUOUS) {
    const [lo, hi] = continuousRange();
    const v = state.cells[state.colorBy][i];
    if (v < 0) return OTHER;
    const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo || 1)));
    return RAMP[Math.round(t * (RAMP.length - 1))];
  }
  const k = categoryOrder(cond).indexOf(categoryOf(i));
  return k >= 0 && k < SLOTS.length ? SLOTS[k] : OTHER;
}

let rangeCache = {};
function continuousRange() {
  const key = state.colorBy;
  if (!rangeCache[key]) {
    const v = state.cells[key].filter(x => x >= 0).sort((a, b) => a - b);
    // 1st-99th percentile so a few extreme cells do not flatten the ramp
    rangeCache[key] = [v[Math.floor(v.length * 0.01)], v[Math.floor(v.length * 0.99)]];
  }
  return rangeCache[key];
}

// ---------- UMAP panels ----------

function buildPanels() {
  const c = state.cells;
  const xs = c.x, ys = c.y;
  const pad = 0.5;
  const xr = [Math.min(...xs) - pad, Math.max(...xs) + pad];
  const yr = [Math.min(...ys) - pad, Math.max(...ys) + pad];
  const host = document.getElementById("panels");

  c.conditions.forEach((cond, k) => {
    const idx = [];
    for (let i = 0; i < state.n; i++) if (c.cond[i] === k) idx.push(i);
    const wrap = document.createElement("div");
    wrap.className = "panel";
    wrap.innerHTML = `<h2>${cond} <span class="n">· ${idx.length.toLocaleString()} cells</span></h2><div class="plot"></div>`;
    host.appendChild(wrap);
    const div = wrap.querySelector(".plot");
    plots.push({ div, idx });

    const axis = { visible: false, fixedrange: false };
    Plotly.newPlot(div, [
      { type: "scattergl", mode: "markers", x: [], y: [], customdata: [], hoverinfo: "none",
        marker: { size: 3, color: [] }, selected: { marker: { opacity: 1 } }, unselected: { marker: { opacity: 1 } } },
      { type: "scattergl", mode: "markers", x: [], y: [], customdata: [], hoverinfo: "none",
        marker: { size: 4.5, color: [] }, selected: { marker: { opacity: 1 } }, unselected: { marker: { opacity: 1 } } },
      // Ring for the cell of the last clicked read. WebGL, because plotly
      // draws WebGL traces above SVG ones.
      { type: "scattergl", mode: "markers", x: [], y: [], hoverinfo: "skip",
        marker: { size: 16, symbol: "circle-open", color: "#0b0b0b", line: { width: 2.5 } } },
    ], {
      margin: { l: 0, r: 0, t: 0, b: 0 }, showlegend: false, dragmode: "lasso",
      xaxis: { ...axis, range: xr }, yaxis: { ...axis, range: yr, scaleanchor: "x" },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      hovermode: "closest",
    }, {
      displaylogo: false, responsive: true,
      modeBarButtonsToRemove: ["toImage", "autoScale2d", "zoomIn2d", "zoomOut2d"],
    });

    div.on("plotly_selected", ev => {
      if (!ev || !ev.points) return;
      applySelection(ev.points.map(p => p.customdata));
      Plotly.relayout(div, { selections: [] });
    });
    div.on("plotly_click", ev => {
      const p = ev.points && ev.points[0];
      if (!p || p.customdata === undefined) return;
      const shift = ev.event && ev.event.shiftKey;
      applySelection([p.customdata], shift ? "add" : undefined);
    });
  });
}

function redrawPanels() {
  const c = state.cells;
  const has = state.selected.size > 0;
  plots.forEach(({ div, idx }, cond) => {
    const bg = { x: [], y: [], cd: [], col: [] };
    const fg = { x: [], y: [], cd: [], col: [] };
    for (const i of idx) {
      const on = !has || state.selected.has(i);
      const t = on && has ? fg : bg;
      t.x.push(c.x[i]); t.y.push(c.y[i]); t.cd.push(i);
      t.col.push(on ? colorFor(i, cond) : DIM);
    }
    const h = state.highlight !== null && c.cond[state.highlight] === cond ? state.highlight : null;
    Plotly.restyle(div, {
      x: [bg.x, fg.x, h === null ? [] : [c.x[h]]],
      y: [bg.y, fg.y, h === null ? [] : [c.y[h]]],
      customdata: [bg.cd, fg.cd, []],
      "marker.color": [bg.col, fg.col, "#0b0b0b"],
      selectedpoints: [null, null, null],
    });
  });
}

// ---------- legend (doubles as group selector) ----------

function renderLegend() {
  const host = document.getElementById("legend");
  host.innerHTML = "";
  const c = state.cells;
  plots.forEach(({ idx }, cond) => {
    const col = document.createElement("div");
    col.className = "col";
    col.innerHTML = `<h3>${c.conditions[cond]}</h3>`;
    if (state.colorBy in CONTINUOUS) {
      const [lo, hi] = continuousRange();
      const fmt = v => (state.colorBy === "umis" ? Math.round(v).toLocaleString() : v.toFixed(2));
      col.innerHTML += `<div class="ramp" style="background:linear-gradient(90deg,${RAMP.join(",")})"></div>
        <div class="ramp-labels"><span>${fmt(lo)}</span><span>${CONTINUOUS[state.colorBy]}</span><span>${fmt(hi)}</span></div>`;
    } else {
      const counts = new Map();
      for (const i of idx) counts.set(categoryOf(i), (counts.get(categoryOf(i)) || 0) + 1);
      const order = categoryOrder(cond).filter(k => counts.has(k));
      for (const k of counts.keys()) if (!order.includes(k)) order.push(k);
      order.forEach((k, slot) => {
        const b = document.createElement("button");
        b.className = "chip";
        const sw = slot < SLOTS.length && categoryOrder(cond).includes(k) ? SLOTS[categoryOrder(cond).indexOf(k)] : OTHER;
        const name = state.colorBy === "cluster" ? `Cluster ${k}` : (k || "unassigned");
        b.innerHTML = `<span class="sw" style="background:${sw}"></span>${name}<span class="n">${counts.get(k).toLocaleString()}</span>`;
        b.title = `Select ${name} in ${c.conditions[cond]}`;
        b.onclick = e => applySelection(idx.filter(i => categoryOf(i) === k), e.shiftKey ? "add" : undefined);
        col.appendChild(b);
      });
    }
    host.appendChild(col);
  });
}

// ---------- selection ----------

function applySelection(indices, modeOverride) {
  const mode = modeOverride || (state.shift ? "add" : state.mode);
  if (mode === "replace") state.selected = new Set(indices);
  else if (mode === "add") indices.forEach(i => state.selected.add(i));
  else indices.forEach(i => state.selected.delete(i));
  onSelectionChanged();
}

function selectedByCondition() {
  const c = state.cells;
  const per = c.conditions.map(() => new Set());
  for (const i of state.selected) per[c.cond[i]].add(c.cb[i]);
  return per;
}

function onSelectionChanged() {
  const c = state.cells;
  const per = selectedByCondition();
  const summary = document.getElementById("sel-summary");
  summary.textContent = state.selected.size === 0
    ? "No selection: reads from all cells"
    : `${state.selected.size.toLocaleString()} cells selected (` +
      c.conditions.map((cond, k) => `${cond} ${per[k].size.toLocaleString()}`).join(" · ") + ")";
  redrawPanels();
  scheduleTrackReload();
}

function exportBarcodes() {
  const c = state.cells;
  const rows = ["condition\tbarcode\tcluster"];
  const ids = state.selected.size ? [...state.selected] : [];
  ids.sort((a, b) => a - b).forEach(i => rows.push(`${c.conditions[c.cond[i]]}\t${c.cb[i]}\t${c.cluster[i]}`));
  const blob = new Blob([rows.join("\n") + "\n"], { type: "text/tab-separated-values" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `${LINE}_selection_${ids.length}cells.tsv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- igv ----------

function readTrackConfig(cond, k, cbs) {
  const c = state.cells;
  const all = state.selected.size === 0;
  const label = all ? "all cells" : `${cbs.size.toLocaleString()} selected cells`;
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
    // downloading the rest of the BAM. Tracks are rebuilt per selection, so
    // the cache buys nothing here anyway.
    cacheBlocks: false,
    displayMode: "SQUISHED", squishedRowHeight: 3,
    showSoftClips: false, showInsertions: false,
    height: 320, coverageTrackHeight: 50,
    order: 10 + k,
  };
  if (state.groupBy !== "none") cfg.groupBy = state.groupBy;
  if (state.readColor !== "none") cfg.colorBy = state.readColor;
  return cfg;
}

let reloadTimer;
function scheduleTrackReload() {
  clearTimeout(reloadTimer);
  reloadTimer = setTimeout(reloadReadTracks, 150);
}

let generation = 0;
async function reloadReadTracks() {
  if (!browser) return; // igv not ready yet; main() loads the tracks once it is
  const gen = ++generation;
  for (const t of readTracks) browser.removeTrack(t);
  readTracks = [];
  const per = selectedByCondition();
  const jobs = state.cells.conditions.map(async (cond, k) => {
    if (state.selected.size && per[k].size === 0) return; // nothing selected here
    try {
      const track = await browser.loadTrack(readTrackConfig(cond, k, per[k]));
      // A newer selection arrived while this track was loading: drop it.
      if (gen !== generation) browser.removeTrack(track);
      else readTracks.push(track);
    } catch (err) {
      console.warn(`Read track for ${cond} failed to load`, err);
    }
  });
  await Promise.all(jobs);
}

async function initIGV() {
  browser = await igv.createBrowser(document.getElementById("igv"), {
    reference: {
      id: "tgme49_toxodb13", name: "T. gondii ME49 (ToxoDB-13)",
      fastaURL: REF + "genome.fasta", indexURL: REF + "genome.fasta.fai",
    },
    locus: START_LOCUS,
    showNavigation: true, showSampleNames: false,
    tracks: [{
      name: "Genes (ToxoDB-13, Cell Ranger GTF)", type: "annotation", format: "gtf",
      url: REF + "genes.gtf.gz", indexURL: REF + "genes.gtf.gz.tbi",
      displayMode: "EXPANDED", order: 1, height: 80,
    }],
  });

  // Clicking a read: find its cell and mark it on the UMAP.
  browser.on("trackclick", (track, data) => {
    if (track.type !== "alignment" || !Array.isArray(data)) return true;
    const cb = data.find(d => d && d.name === "CB");
    if (!cb) return true;
    const k = readTracks.indexOf(track) >= 0 ? state.cells.conditions.indexOf(
      Object.keys(state.cells.bams).find(cond => track.config.url.endsWith(state.cells.bams[cond]))) : -1;
    const i = k >= 0 ? cbIndex[k].get(String(cb.value)) : undefined;
    if (i !== undefined) showReadCell(i);
    return true;
  });
}

function showReadCell(i) {
  const c = state.cells;
  state.highlight = i;
  const host = document.getElementById("read-cell");
  host.innerHTML = `Read from <b>${c.cb[i]}</b> · ${c.conditions[c.cond[i]]} cluster ${c.cluster[i]}
    · ${c.umis[i].toLocaleString()} UMIs <button>Select cell</button>`;
  host.querySelector("button").onclick = () => applySelection([i], "replace");
  redrawPanels();
}

// ---------- gene search ----------

let genes = [];
function initSearch() {
  const input = document.getElementById("search-input");
  const list = document.getElementById("search-results");
  const go = locus => { list.hidden = true; browser.search(locus); };
  const locusOf = g => {
    const pad = Math.round((g[3] - g[2]) * 0.15) + 200;
    return `${g[1]}:${Math.max(1, g[2] - pad)}-${g[3] + pad}`;
  };
  const matches = q => {
    q = q.trim().toLowerCase();
    if (q.length < 2) return [];
    const exact = genes.filter(g => g[0].toLowerCase() === q);
    const rest = genes.filter(g => g[0].toLowerCase() !== q &&
      (g[0].toLowerCase().includes(q) || g[5].toLowerCase().includes(q)));
    return exact.concat(rest).slice(0, 50);
  };
  input.addEventListener("input", () => {
    const m = matches(input.value);
    list.innerHTML = m.map((g, j) =>
      `<li data-j="${j}"><span class="gid">${g[0]}</span><span class="desc">${g[5]}</span></li>`).join("");
    list.hidden = m.length === 0;
    list.querySelectorAll("li").forEach(li => li.onclick = () => {
      const g = m[+li.dataset.j]; input.value = g[0]; go(locusOf(g));
    });
  });
  document.getElementById("search").addEventListener("submit", e => {
    e.preventDefault();
    const m = matches(input.value);
    go(m.length ? locusOf(m[0]) : input.value.trim());
  });
  document.addEventListener("click", e => { if (!e.target.closest("#search")) list.hidden = true; });
}

// ---------- wiring ----------

async function main() {
  const [cells, geneList] = await Promise.all([loadJSON(DATA + "cells.json"), loadJSON(REF + "genes.json")]);
  state.cells = cells;
  state.n = cells.cb.length;
  genes = geneList;
  cbIndex = cells.conditions.map(() => new Map());
  for (let i = 0; i < state.n; i++) cbIndex[cells.cond[i]].set(cells.cb[i], i);
  document.getElementById("line-label").textContent = cells.line;

  buildPanels();
  redrawPanels();
  renderLegend();
  await initIGV();
  initSearch();
  await reloadReadTracks();

  document.getElementById("color-by").onchange = e => {
    state.colorBy = e.target.value; redrawPanels(); renderLegend();
  };
  document.getElementById("group-by").onchange = e => { state.groupBy = e.target.value; scheduleTrackReload(); };
  document.getElementById("read-color").onchange = e => { state.readColor = e.target.value; scheduleTrackReload(); };
  document.querySelectorAll(".seg button").forEach(b => b.onclick = () => {
    state.mode = b.dataset.mode;
    document.querySelectorAll(".seg button").forEach(x => x.classList.toggle("on", x === b));
  });
  document.getElementById("sel-clear").onclick = () => { state.selected = new Set(); onSelectionChanged(); };
  document.getElementById("sel-export").onclick = exportBarcodes;
  addEventListener("keydown", e => { if (e.key === "Shift") state.shift = true; });
  addEventListener("keyup", e => { if (e.key === "Shift") state.shift = false; });
}

main().catch(err => {
  document.body.insertAdjacentHTML("afterbegin",
    `<p style="padding:12px;background:#fde8e8;margin:0">Failed to load: ${err.message}.
     Serve the project with <code>python3 scripts/serve.py</code> and run the prep scripts first.</p>`);
  throw err;
});
