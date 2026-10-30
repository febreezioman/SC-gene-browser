# Session — 2026-10-30: project setup and Phase 1 prototype

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

## Follow-up in the same session: Phase 1 built
- User answers: IAA timing unknown to them → found in Sc-x2 overview (24h);
  prototype on AP2XI-2/AP2XII-1; local first, public later; include antisense;
  single-cell only.
- Wrote `scripts/config.py`, `prep_ref.py`, `prep_bam.py`, `export_cells.py`,
  `serve.py`, and `web/` (index.html, app.js, style.css).
- Found igv.js accepts a function filter → no fork needed.
- Found and worked around an igv.js 3.8.8 block-cache bug (open-ended Range).
- Fixed a track-reload race (generation counter, parallel loads) and the read
  highlight ring (must be WebGL to sit above WebGL points).
- Verified UMI counts match Cell Ranger exactly at two genes.

## Follow-up 2: usability round (user feedback)
- Selection is staged: build it (click toggles a whole group, legend chips
  toggle, lasso adds/removes), then **Load reads**. No reads load while selecting.
  Chips show full/partial selection; the summary lists the groups.
- Reads: default no grouping (one pile per condition), expanded rows, tracks
  fit the panel height; Track height slider + Fit; Display expanded/squished.
- Draggable splitter between the panels (width remembered in localStorage).
- Gene models labelled by ME49 gene ID (transcript_id rewritten; genes and
  transcripts are 1:1). Search by common names pulled from ToxoDB-13/68
  descriptions (1,593 genes), ignoring case and hyphens.
- igv.js has no public resize(); the page calls its window-resize handler.

## Follow-up 3: GitHub, Crucial X10, AlRubaye line
- Pushed to github.com/febreezioman/SC-gene-browser (private).
- Copied data/ to the Crucial X10 (exFAT: no symlinks, so FASTAs are copied).
- Generalised the pipeline: config REFS + LINES (kind cellranger/starsolo),
  fields declared per line in cells.json, line switcher, data/ref/<ref>/.
- Added AlRubaye 2026 cat gut (SRR36982562, STARsolo, ToxoDB-68) with our own
  UMI dedup; verified against STARsolo counts.
- Slow writes: 16 parallel writers to the NTFS data drive stalled (~1 MB/s
  each); per-contig parts now go to ext4 scratch (229 MB/s). The final write
  then queued behind a large Thunar copy to the same drive (41 GB dirty pages).

## Decisions made
- igv.js 3.8.8 without a fork; Sc-x2 read-only; generated data in `data/`.
- Default locus AP2XI-2 (SAG1 is too dense to start on).

## Next session
- Try it by hand; collect UX wishes.
- Decide on the Phase 2 backend (also the route to public hosting).
- Add more lines.
