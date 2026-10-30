"""Paths and the dataset -> BAM mapping. Sc-x2 is read-only input."""
from pathlib import Path

SCX2 = Path.home() / "projects" / "Sc-x2"
CR = SCX2 / "raw" / "CellRanger"
REF = SCX2 / "refs" / "Reference-Genome-Toxo13"
COORDS = SCX2 / "results" / "figs" / "umap3d" / "umap3d_coords.csv.gz"

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

GENOME_FASTA = REF / "ToxoDB-13.0_TgondiiME49_Genome.fasta"
GTF = REF / "Tgondii_ToxoDB13.cellranger.gtf"
ANNOTATION_TSV = REF / "ToxoDB-68_ME49_master_annotation.tsv"

# One entry per line. `samples` maps a condition (as written in COORDS) to its
# Cell Ranger run. Labels use AP2 gene names (see AI.md).
LINES = {
    "AP2XI-2_AP2XII-1": {
        "label": "AP2XI-2/AP2XII-1",
        "coords_dataset": "XI-2_XII-1",
        "h5ad": SCX2 / "results" / "figs" / "duo_AP2XI-2_AP2XII-1" / "combined.h5ad",
        "samples": {
            "UT": CR / "Duo_UT_Feb2026",
            # Folder says 48h but the run is IAA 24h: relabelled in Sc-x2 on
            # 2026-06-30, folder name kept (Sc-x2 docs/overview.md).
            "IAA 24h": CR / "Duo_IAA_48h_Feb2026",
        },
    },
}


def sample_id(condition: str) -> str:
    """File-safe id for a condition, e.g. 'IAA 24h' -> 'IAA_24h'."""
    return condition.replace(" ", "_")


def bam_cb(barcode: str) -> str:
    """Table barcode 'AAAC...-1-0' -> BAM CB tag 'AAAC...-1'."""
    return barcode.rsplit("-", 1)[0]
