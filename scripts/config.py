"""Paths, references and the line -> BAM mapping. All inputs are read-only."""
from pathlib import Path

SCX2 = Path.home() / "projects" / "Sc-x2"
CR = SCX2 / "raw" / "CellRanger"
COORDS = SCX2 / "results" / "figs" / "umap3d" / "umap3d_coords.csv.gz"
GENOME_FILES = Path.home() / "projects" / "Genome-files"
EES_REFS = Path("/media/jon/Working-drive-nv/EES-Oocyst-RNAseq/refs/ME49")
ALRUBAYE = SCX2 / "processing" / "alrubaye2026_perturbseq"
ALRUBAYE_CELLS = SCX2 / "results" / "figs" / "novel_antisense_10x" / "alrubaye_perturbseq"

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
# Local ext4 scratch for prep_bam.py's per-contig parts: many parallel writers
# to the NTFS drive that holds data/ stall (~1 MB/s each).
SCRATCH = Path("/media/jon/Working-drive-nv/scratch/sc-read-browser")

_TOXO13 = SCX2 / "refs" / "Reference-Genome-Toxo13"
# Gene descriptions and common names (same gene IDs in both releases).
ANNOTATION_TSV = _TOXO13 / "ToxoDB-68_ME49_master_annotation.tsv"
NAMES_GFF = _TOXO13 / "ToxoDB-13.0_TgondiiME49.gff"

# Each reference must match the contigs of the BAMs aligned to it.
REFS = {
    "ToxoDB-13": {
        "label": "T. gondii ME49 (ToxoDB-13)",
        "fasta": _TOXO13 / "ToxoDB-13.0_TgondiiME49_Genome.fasta",
        "gtf": _TOXO13 / "Tgondii_ToxoDB13.cellranger.gtf",  # what Cell Ranger counted against
    },
    "ToxoDB-68": {
        "label": "T. gondii ME49 (ToxoDB-68)",
        "fasta": GENOME_FILES / "ToxoDB-68" / "ToxoDB-68_TgondiiME49_Genome.fasta",
        "gtf": EES_REFS / "ToxoDB-68_TgondiiME49.gtf",  # in the STARsolo index
    },
}

# One entry per line (dataset). `samples` maps a condition to its BAM.
# kind "cellranger": duplicates flagged 0x400, CB has a "-1" suffix.
# kind "starsolo": no duplicate flags (deduplicated by prep_bam.py), CB bare.
# Labels use AP2 gene names (see AI.md).
LINES = {
    "AP2XI-2_AP2XII-1": {
        "label": "AP2XI-2/AP2XII-1",
        "kind": "cellranger",
        "ref": "ToxoDB-13",
        "coords_dataset": "XI-2_XII-1",
        "h5ad": SCX2 / "results" / "figs" / "duo_AP2XI-2_AP2XII-1" / "combined.h5ad",
        "samples": {
            "UT": CR / "Duo_UT_Feb2026" / "outs" / "possorted_genome_bam.bam",
            # Folder says 48h but the run is IAA 24h: relabelled in Sc-x2 on
            # 2026-06-30, folder name kept (Sc-x2 docs/overview.md).
            "IAA 24h": CR / "Duo_IAA_48h_Feb2026" / "outs" / "possorted_genome_bam.bam",
        },
    },
    "AlRubaye2026_cat_gut": {
        "label": "AlRubaye 2026 cat gut (enteric stages)",
        "kind": "starsolo",
        "ref": "ToxoDB-68",
        # SRR36982562, STARsolo on adapter-trimmed reads against ME49 + cat +
        # human + mouse (Sc-x2 scripts/58); this BAM is the unique ME49 reads.
        "samples": {"Cat gut": ALRUBAYE / "ME49_unique.bam"},
    },
}


def sample_id(condition: str) -> str:
    """File-safe id for a condition, e.g. 'IAA 24h' -> 'IAA_24h'."""
    return condition.replace(" ", "_")


def bam_cb(barcode: str) -> str:
    """Cell Ranger table barcode 'AAAC...-1-0' -> BAM CB tag 'AAAC...-1'."""
    return barcode.rsplit("-", 1)[0]
