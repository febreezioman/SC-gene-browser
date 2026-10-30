# Data inventory

All inputs are **read-only**. Paths are under `~/projects/Sc-x2` unless stated.
Surveyed 2026-10-30.

## Alignments: Cell Ranger BAMs

`raw/CellRanger/<sample>/outs/possorted_genome_bam.bam` (+ `.bai`).
Cell Ranger 10.0.0, coordinate sorted, 11–14 GB each.

| Sample folder | Line | Condition (coords table) |
|---|---|---|
| `SC_A2_AP2XII1_KD_AP2VIII7_KD_UT` | AP2XII-1/AP2VIII-7 | UT |
| `SC_B2_AP2XII1_KD_AP2VIII7_KD_IAA-48h` | AP2XII-1/AP2VIII-7 | IAA 48h |
| `Duo_UT_Feb2026` | AP2XI-2/AP2XII-1 | UT |
| `Duo_IAA_48h_Feb2026` | AP2XI-2/AP2XII-1 | **IAA 24h** (folder name is wrong, see below) |
| `SC_E2_TripleKD_UT` | AP2XII-1/AP2VIII-7/AP2XI-2 | UT |
| `SC_F2_TripleKD_IAA_24h` | AP2XII-1/AP2VIII-7/AP2XI-2 | IAA 24h |
| `SC_G2_TripleKD_IAA_48h` | AP2XII-1/AP2VIII-7/AP2XI-2 | IAA 48h |
| `SC_C2_MORC-KD_UT` | MORC KD | UT (not in the coords table) |
| `SC_D2_MORC-KD_IAA-24h` | MORC KD | IAA 24h (not in the coords table) |

**AP2XI-2/AP2XII-1 IAA is 24h.** The folder name says 48h because of a
mislabel at data organisation. Sc-x2 corrected the label on 2026-06-30 and
kept the folder name (Sc-x2 `docs/overview.md`, 2026-06-30 entry). Still stale
in Sc-x2: the sample table further down that same overview, and the `label`
column of `duo_AP2XI-2_AP2XII-1/combined.h5ad` (`IAA-48h`). This project uses
its own labels from `scripts/config.py`.

Read tags present (from a `Duo_UT_Feb2026` sample):
`CB` corrected barcode (`TACACCTTCGACCGCA-1`), `UB` corrected UMI,
`CR`/`UR` raw, `xf` (count-flag bitfield; reads used for counting carry
`xf & 8`), `AN` antisense alignment, `RE` region (E/N/I), `NH`/`HI` multimap,
`RG` = one value per library (not per cell). Duplicate reads are flagged 0x400.

## Cells, embeddings, clusters

- `results/figs/umap3d/umap3d_coords.csv.gz`: **best single source.** One row
  per cell: `dataset, condition, barcode, UMAP3D_1..3, leiden_cluster,
  phase_call, pseudotime, AUCell_Bradyzoite, AUCell_EES1to3, brady_ees_call,
  ees_10stage_call`. 7 dataset×condition groups, 112,685 cells.
- `results/figs/<line>/combined.h5ad`: per-line AnnData with `X_umap` (2D),
  `X_pca`, `leiden`, `phase_call`, cell-cycle scores. Needed for **2D** UMAP
  (lasso selection) and for gene expression colouring.
- `drafts/Single-Cell-Draft-Figs/clustering/leiden_cell_assignments.csv`:
  locked resolution-0.4 Leiden clusters (per condition, so cluster 0 in UT ≠
  cluster 0 in IAA).
- Existing page to build from: `results/figs/umap3d/umap3d_interactive.html`
  (from `scripts/52_umap3d_interactive.py`, plotly 3D, ~10 MB).

**Barcode format:** tables use `AAACCAAAGATAGCCC-1-0` (the trailing `-0` is the
AnnData concat index). The BAM `CB` is `AAACCAAAGATAGCCC-1`. Strip the last
`-N` and pair with the right BAM via dataset + condition.

## Reference

`refs/Reference-Genome-Toxo13/`:
- `ToxoDB-13.0_TgondiiME49_Genome.fasta`: 14 `TGME49_chr*` chromosomes plus
  ~2,250 small unplaced contigs (`tgme49_asmbl.*`, `ABPA02*`), ~65 Mb. 98% of
  reads land on the chromosomes. Gene search is still the easiest way around.
- `Tgondii_ToxoDB13.cellranger.gtf`: what the BAMs were counted against
  (8,920 genes; transcript/exon/CDS, `gene_id` only). Used as the annotation
  track. Contig names match the BAM header.
- `ToxoDB-68_ME49_master_annotation.tsv`: gene descriptions for search.
  SAG1 is listed as "SAG-related sequence SRS29B" (TGME49_233460).

## Other data

- `processing/alrubaye2026_perturbseq/`: STARsolo BAM from a public dataset,
  plus per-cell-type BAMs (`celltype/*.bam`). Possible second dataset later.
- `/media/jon/Working-2-SSD/EES-Oocyst-RNAseq_tmp/coug_flair2/`: bulk
  long-read (FLAIR) oocyst data. **Out of scope** (user, 2026-10-30): the project
  is single-cell only for now.

## Generated files (`data/`, gitignored)

| File | Made by | Notes |
|---|---|---|
| `data/ref/genome.fasta(.fai)` | `prep_ref.py` | symlink to Sc-x2 FASTA + index |
| `data/ref/genes.gtf.gz(.tbi)` | `prep_ref.py` | sorted, bgzipped, `gene_name` = gene ID |
| `data/ref/genes.json` | `prep_ref.py` | gene search index with descriptions |
| `data/<line>/<condition>.bam(.bai)` | `prep_bam.py` | called cells only, no 0x400 dups, `CL` tag |
| `data/<line>/cells.json` | `export_cells.py` | 2D UMAP + per-cell metadata, columnar |
