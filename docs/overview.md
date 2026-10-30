# Overview

## Goal

Select cells on an embedding (one cell, lasso, cluster, metadata filter) and see
the gene annotation and the exact reads those cells contributed, in a genome
browser view. Prototype first as a "souped-up HTML file", then decide whether
it needs a backend or a fork of igv.js.

## Status (2026-10-30)

Two datasets ("lines") in one page, switched from the header:
- AP2XI-2/AP2XII-1 (Cell Ranger, ToxoDB-13): UT 11,642 + IAA 24h 14,386 cells.
- AlRubaye 2026 cat gut (STARsolo, ToxoDB-68): 23,121 enteric-stage cells,
  95.5M reads after dedup (340M duplicate reads dropped, ~4.6 reads/molecule).
  Dedup checked against STARsolo counts: 495/495, 756/756, 2,316/2,320 UMIs.

Code is on GitHub (private): febreezioman/SC-gene-browser. A copy of `data/`
is on the Crucial X10 (`SC-gene-browser/data/`) for use on another computer.

- Slim BAMs: UT 28.2M reads / 1.50 GB, IAA 24h 49.4M reads / 2.71 GB
  (from 11–12 GB each), ~7 min per condition.
- Verified: the slim BAM holds exactly the molecules Cell Ranger counted
  (AP2XI-2: 222/222 UMIs, AP2XII-1: 985/985 in UT).
- Verified in headless Chrome: cluster selection, shift-add across conditions,
  gene search, read click → cell ring. AP2XI-2 locus loads ~6 MB per track.
- Known limit: every read in the window is downloaded and parsed before the
  cell filter. Fine for most genes; slow for very high ones (SRS29B/SAG1 has
  ~583k UT reads in 2 kb). That is the trigger for the Phase 2 backend.
- Observation: most UT reads at AP2XI-2 are antisense to it (`AN` tag, + strand
  on a − strand gene), and Cell Ranger did not count them (`xf 0`).

## Decisions

- **igv.js 3.8.8, no fork.** Its alignment `filter` accepts a function, so the
  page filters on the `CB` tag directly. No read-group rewrite.
- `cacheBlocks: false` on read tracks to avoid an igv.js 3.8.8 bug (open-ended
  Range request that pulled 790 MB per load). Worth reporting upstream.
- Slim BAMs keep everything from called cells except UMI duplicates, including
  antisense, intronic, intergenic and multimapped reads (user: antisense must be
  visible).
- AP2XI-2/AP2XII-1 IAA is **24h**. The `Duo_IAA_48h_Feb2026` folder name is a
  known mislabel (Sc-x2 overview, 2026-06-30).
- Single-cell only for now; the EES oocyst long reads are out of scope.
- Local use first; public hosting is wanted later (see design.md, Hosting).
- Sc-x2 is read-only; all generated data lives in `data/` (gitignored).
- 2D UMAP for selection, one panel per condition (clusters are per condition).
- Selection is staged and reads load on demand (user: live reloading was laggy).
- One read track per condition, reads not grouped by default. Merging the
  conditions into one track would need a merged BAM with a condition tag
  (barcodes can repeat across conditions); not done.

- Each line's cells.json declares its fields (colour/group options), so new
  datasets need no page changes; prep_bam.py takes barcodes and CL tags from it.
- STARsolo dedup: one read per (CB, UB, GX gene), else per (CB, UB, strand,
  start). Parallel per contig, parts on ext4 scratch (NTFS stalls).
- Gene models labelled via igv.js `nameField: "gene_id"` (keeps isoforms).

## Open questions / next

1. Phase 2 backend (FastAPI + pysam returning only the selected cells' reads),
   so highly expressed genes and large selections stay fast. Needed for public
   hosting too, since that also needs real auth.
2. Add the other AP2 lines (AP2XII-1/AP2VIII-7, AP2XII-1/AP2VIII-7/AP2XI-2, MORC KD).
   MORC KD has no rows in the coords table; its cells/clusters need another source.
3. Per-cell normalised coverage: selection vs rest.
4. Colour the UMAP by expression of the gene in view.
