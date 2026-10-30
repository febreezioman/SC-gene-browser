# Session — 2026-10-30: project setup

## Goal
Set up the project and docs for a single-cell browser that shows per-cell reads.

## What was done
- Created README, AI.md/CLAUDE.md protocol (modelled on Sc-x2), docs, layout.
- Surveyed inputs in `~/projects/Sc-x2` and the EES oocyst folder → `docs/data.md`.
- Wrote the architecture and phased plan → `docs/design.md`.

## Results / findings
- BAMs carry `CB`/`UB`/`xf`/`AN`/`RE` tags; contigs match the Cell Ranger GTF.
- Table barcodes have a `-N` sample suffix that the BAM `CB` lacks.
- `umap3d_coords.csv.gz` has every cell with clusters and stage calls (3D only;
  2D UMAP comes from the per-line `combined.h5ad`).
- EES oocyst data is bulk long reads (FLAIR), not single-cell.
- Possible label mismatch: AP2XI-2/AP2XII-1 IAA is 48h (folder) vs 24h (table).

## Decisions made
- igv.js, no fork for Phase 1; Sc-x2 read-only; generated data in `data/`.

## Next session
- Answer the open questions in `docs/overview.md`.
- Write `scripts/prep_bam.py` on one dataset; measure the filtered BAM size.
- Write `scripts/export_cells.py`, then the first `web/index.html`.
