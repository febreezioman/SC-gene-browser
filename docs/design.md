# Design

## The idea

Embedding on the left, genome browser on the right, one linked state:

```
┌──────────────────────┐ ┌──────────────────────────────────────┐
│  UMAP (2D, plotly)   │ │ locus: gene search / coords          │
│  colour: cluster ▾   │ │ ── gene models (GTF) ─────────────── │
│  click / lasso /     │ │ ── coverage: selection vs. rest ──── │
│  legend = selection  │ │ ── reads: only selected cells,       │
│                      │ │    grouped/coloured by cluster ───── │
└──────────────────────┘ └──────────────────────────────────────┘
   selection: 214 cells (UT · 3, lasso)  [clear] [export barcodes]
```

Selection sources: single click, box/lasso (2D only; plotly cannot lasso in
3D), legend/cluster picker, or a metadata filter (phase, stage call, ...).

## Tool choice

| Option | License | Read-level? | Filter by CB list | Verdict |
|---|---|---|---|---|
| **igv.js** | MIT | yes | via RG trick, or small fork | **Start here** |
| JBrowse 2 | Apache-2.0 | yes | `tags.CB` field filter (v5 beta, untested) | Revisit if igv.js limits bite |
| JBrowse single-cell demo | Apache-2.0 | coverage only | whole cell types | Too coarse |
| Vitessce | MIT | no | n/a | Profiles only |
| UCSC / IGB | mixed / Java desktop | n/a | n/a | Skip |

igv.js is one JS file, no build step, and the CB filter is a small,
self-contained change (and upstreamable).

## Phase 1: static prototype, no fork

1. **`scripts/prep_bam.py`** (pysam). Per dataset×condition:
   - keep reads whose `CB` is in the cell table, drop 0x400 duplicates
     (one read per molecule), optionally keep only `xf & 8` (counted reads);
   - **rewrite `RG` = `<dataset>|<condition>|<CB>`** so each cell is a read
     group, and add a cluster tag (e.g. `CL:Z:UT·3`);
   - write `data/<id>.bam` + index. Expect a large shrink from 11–14 GB;
     measure it.
2. **`scripts/export_cells.py`**: one compact cell table for the page
   (barcode as in the BAM RG, 2D UMAP, cluster, phase, stage calls).
   JSON or Arrow; ~110k cells is fine in the browser.
3. **`web/index.html`**: plotly 2D UMAP + igv.js. On selection, set the
   alignment track's `filter.readgroups` to the selected cells, `groupBy:
   tag:CL`, raise `samplingDepth` well above the default 100.
4. **Local server with Range support** (e.g. `npx http-server`, or a 20-line
   Python server). Plain `python -m http.server` will not work.

Done when: clicking one cell at a highly expressed gene (e.g. SAG1) shows its
reads, and the count matches that cell's UMI count in the h5ad.

## Phase 2: things Phase 1 will expose

- **Everything in the window is downloaded, then filtered.** For a whole
  dataset that may be fine at one gene on a 65 Mb genome, but measure.
  If not: a tiny local backend (FastAPI + pysam) that returns only the selected
  cells' reads for a region, as a BAM/SAM slice or JSON. Scales to any selection
  and also enables the coverage comparison below cheaply.
- **Selection vs. rest coverage**, normalised per cell. Compute server-side
  or in the page from the loaded reads.
- **Fork igv.js** for a real `filter: {tag: 'CB', values: [...]}` instead of
  the RG rewrite; offer it upstream.
- Gene expression colouring of the UMAP, linked to the current locus.
- Multiple datasets side by side (UT vs IAA tracks for the same selection).

## Hosting (later)

- GitLab Pages: range requests only with `ARTIFACT_COMPRESSION_LEVEL:
  fastest`; 1 GB site cap, so the full BAMs will not fit. A per-gene subset
  might.
- StatiCrypt only encrypts the HTML. BAMs next to it are public to anyone with
  the URL. Unpublished data needs real auth (basic-auth server) or local-only use.

## Biology-specific things to remember

- 10x 3' chemistry: reads pile up at 3' ends; transcript bodies are sparse.
  Show strand, and expect antisense (`AN` tag) and intronic (`RE:I`) reads,
  which Sc-x2 already studies (novel antisense scripts 53–61).
- Fragmented assembly (2,263 contigs): gene-ID search is the main navigation.
