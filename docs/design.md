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
| **igv.js** | MIT | yes | **yes: `filter` accepts a function** (no fork) | **Using it (3.8.8)** |
| JBrowse 2 | Apache-2.0 | yes | `tags.CB` field filter (v5 beta, untested) | Revisit if igv.js limits bite |
| JBrowse single-cell demo | Apache-2.0 | coverage only | whole cell types | Too coarse |
| Vitessce | MIT | no | n/a | Profiles only |
| UCSC / IGB | mixed / Java desktop | n/a | n/a | Skip |

igv.js is one JS file with no build step. **Found 2026-10-30 in the 3.8.8
source:** an alignment track's `filter` may be a function
(`typeof config.filter === 'function'` in `setReaderDefaults`), and it runs
*before* downsampling. So `filter: a => selectedCBs.has(a.getTag('CB'))` works
on the stock library: no read-group rewrite, no fork. A function filter
replaces the default filter, so the page re-applies the default checks itself
(mapped, not QC-fail, not duplicate).

## Phase 1: static prototype, no fork

1. **`scripts/prep_bam.py`**: per condition, `samtools view -D CB:<cells>
   -F 0x404` keeps called cells only and drops UMI duplicates (one read per
   molecule). **All other reads are kept, including antisense (`AN`), intronic
   and intergenic ones and multimappers** (user decision: antisense must be
   visible). pysam then adds `CL:Z:<leiden cluster>` for grouping.
2. **`scripts/prep_ref.py`**: indexed FASTA, bgzipped+tabixed GTF, gene search
   JSON with ToxoDB descriptions.
3. **`scripts/export_cells.py`**: columnar `cells.json` (2D UMAP from
   `combined.h5ad`, clusters/phase/stage calls from the coords table).
4. **`web/`** (`index.html`, `app.js`, `style.css`): one plotly panel per
   condition on the shared 2D embedding (clusters are per condition, so never
   on a shared key). Lasso/click/legend-chip selection with replace/add/remove.
   One igv.js read track per condition, rebuilt on selection change, grouped by
   `tag:CL`, coloured by strand / region (`RE`) / cell (`CB`),
   `samplingDepth: 10000`. Clicking a read rings its cell on the UMAP.
5. **`scripts/serve.py`**: static server with Range support. Plain
   `python -m http.server` ignores Range.

Run:
```bash
P=~/projects/Sc-x2/envs/sc/bin/python
$P scripts/prep_ref.py
$P scripts/prep_bam.py AP2XI-2_AP2XII-1
$P scripts/export_cells.py AP2XI-2_AP2XII-1
python3 scripts/serve.py        # http://127.0.0.1:8000/web/
```

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
- A fork of igv.js is probably unnecessary (see function filter above). Fork
  only for things like per-cell normalised coverage in the track itself.
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
