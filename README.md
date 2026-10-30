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

## Status

Planning. Nothing runs yet; see `docs/design.md` for Phase 1.
