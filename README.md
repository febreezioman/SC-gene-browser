# Single-cell browser with reads

Prototype of a single-cell browser that links an embedding (UMAP/PCA) to a
genome browser. Pick one cell, a lasso of cells, or whole clusters, and see the
gene annotation plus the **actual aligned reads** those cells contributed,
in an IGV/JBrowse-style view.

First target data: the *T. gondii* 10x 3' datasets in `~/projects/Sc-x2`
(AP2 knockdown lines, ToxoDB-13 ME49 reference).

See [`AI.md`](AI.md) for the AI working protocol, and
[`docs/overview.md`](docs/overview.md) for goals, status and decisions.

## Layout

```
docs/
  overview.md      # Living overview: goals, status, decisions, open questions
  design.md        # Architecture, tool choice, phased plan
  data.md          # Input data inventory (BAMs, h5ad, refs) and caveats
  sessions/        # Per-session logs
scripts/           # Python prep: BAM subsetting/retagging, cell table export
web/               # The browser page (HTML/JS), igv.js + plotly
data/              # Generated, browser-ready files (gitignored)
```

## Run it

```bash
P=~/projects/Sc-x2/envs/sc/bin/python
$P scripts/prep_ref.py                          # reference, ~1 min
$P scripts/prep_bam.py AP2XI-2_AP2XII-1         # slim BAMs, ~15 min
$P scripts/export_cells.py AP2XI-2_AP2XII-1     # cell table
python3 scripts/serve.py                        # then open http://127.0.0.1:8000/web/
```

Lasso/click cells or click a legend entry; reads for those cells appear on the
right, one track per condition. Click a read to ring its cell on the UMAP.
`?locus=TGME49_chrXII:1176841-1187769` opens at a given locus.

## Status

Phase 1 prototype works on AP2XI-2/AP2XII-1 (UT, IAA 24h). See `docs/overview.md`.
