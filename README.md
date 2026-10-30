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

Viewing needs only Python 3 and the `data/` folder:

```bash
python3 scripts/serve.py        # then open http://127.0.0.1:8000/web/
```

Switch datasets with the menu in the header (or `?line=<id>`). Click a cell or
a legend entry to toggle its whole group (or lasso cells), then press
**Load reads**: one read track per condition appears on the right. Drag the bar
between the panels to resize; read tracks fit the panel height (or use the
Track height slider). Search by common name (SAG1, AP2XI-2), TGME49 ID or
description. Click a read to ring its cell on the UMAP.
`?locus=TGME49_chrXII:1176841-1187769` opens at a given locus.

Building data (only on the machine with Sc-x2 and the raw BAMs; paths in
`scripts/config.py`):

```bash
P=~/projects/Sc-x2/envs/sc/bin/python
$P scripts/prep_ref.py                      # references, ~1 min
$P scripts/export_cells.py <line>           # cells.json + data/lines.json
$P scripts/prep_bam.py <line>               # slim BAMs (reads cells.json)
```

Lines: `AP2XI-2_AP2XII-1` (Cell Ranger, ToxoDB-13) and `AlRubaye2026_cat_gut`
(STARsolo, ToxoDB-68).

## Status

Working prototype with two datasets. See `docs/overview.md`.
