# Overview

## Goal

Select cells on an embedding (one cell, lasso, cluster, metadata filter) and see
the gene annotation and the exact reads those cells contributed, in a genome
browser view. Prototype first as a "souped-up HTML file", then decide whether
it needs a backend or a fork of igv.js.

## Status (2026-10-30)

- Project set up: docs, protocol, layout. No code yet.
- Data surveyed (`docs/data.md`): 9 Cell Ranger BAMs, a per-cell table with
  112,685 cells across 7 dataset×condition groups, ToxoDB-13 reference.
- Plan written (`docs/design.md`): Phase 1 = igv.js + plotly, BAMs retagged
  so each cell is a read group.

## Decisions

- **igv.js** over JBrowse 2 for the prototype: MIT, single file, easy to fork,
  stable. JBrowse kept as a fallback.
- Start **without a fork** (RG = cell trick); fork only if needed.
- Sc-x2 is read-only; all generated data lives in `data/` (gitignored).
- 2D UMAP for selection (no lasso in plotly 3D). The 3D page can come later.

## Open questions

1. AP2XI-2/AP2XII-1 IAA: the BAM folder says 48h, the coords table says 24h.
   Which is right?
2. First dataset to prototype on? Suggest AP2XI-2/AP2XII-1 UT + IAA (smallest).
3. Local-only use, or does it need to go on the lab site? (Decides hosting/auth.)
4. Count only reads used for UMI counting (`xf & 8`), or show everything,
   including intronic/antisense reads that were not counted?
5. EES oocyst long reads: wanted as a track (isoforms / aligned reads)?
