#!/usr/bin/env python3
"""
Browser-ready reference: genome FASTA + index, annotation GTF, gene search index.

  data/ref/genome.fasta(.fai)  symlink to the Sc-x2 FASTA, indexed here
  data/ref/genes.gtf.gz(.tbi)  the Cell Ranger GTF with transcript_id and
                               gene_name set to the ME49 gene ID, so igv.js
                               labels models "TGME49_233460" (not
                               "rna_TGME49_233460-1"; genes and transcripts are
                               1:1 here). Sorted, bgzipped, tabix-indexed.
  data/ref/genes.json          [[gene_id, chr, start, end, strand, description,
                                 [common names]]] for the page's gene search

Common names are not a field in either ToxoDB file; they are pulled from the
descriptions: parenthesised symbols in ToxoDB-13 ("SRS29B (SAG1)") and
symbol-like trailing words in ToxoDB-68 ("... AP2XI-2"). ~1,600 genes get one.

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/prep_ref.py
"""
from __future__ import annotations

import json
import re
import urllib.parse

import pandas as pd
import pysam

import config

# A symbol has a capital letter and a digit: SAG1, AP2XI-2, ROP40, ApiAT7-1.
SYMBOL = re.compile(r"^(?=.*\d)(?=.*[A-Z])[A-Za-z][A-Za-z0-9]*(?:[-.][A-Za-z0-9]+)*$")


def symbols(desc: str) -> list[str]:
    out = []
    for inner in re.findall(r"\(([^()]+)\)", desc):
        out += [t for t in re.split(r"[,;/ ]+", inner) if SYMBOL.match(t)]
    words = re.sub(r"\([^()]*\)", "", desc).replace(",", " ").split()
    out += [t for t in words[-2:] if SYMBOL.match(t) and len(t) <= 12]
    return out


def common_names() -> dict[str, list[str]]:
    names: dict[str, list[str]] = {}
    with open(config.GFF) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                attrs = dict(a.split("=", 1) for a in f[8].strip().split(";") if "=" in a)
                desc = urllib.parse.unquote_plus(attrs.get("description", ""))
                names.setdefault(attrs["ID"], []).extend(symbols(desc))
    ann = pd.read_csv(config.ANNOTATION_TSV, sep="\t", usecols=["gene_id", "description"])
    for gid, desc in zip(ann.gene_id, ann.description.fillna("")):
        names.setdefault(gid, []).extend(symbols(desc))
    # Keep the hyphenated spelling when both "AP2XI2" and "AP2XI-2" occur.
    clean = {}
    for gid, syms in names.items():
        uniq = list(dict.fromkeys(syms))
        hyphenated = {s.replace("-", "") for s in uniq if "-" in s}
        clean[gid] = [s for s in uniq if "-" in s or s not in hyphenated]
    return clean


def main() -> None:
    out = config.DATA / "ref"
    out.mkdir(parents=True, exist_ok=True)

    fasta = out / "genome.fasta"
    if not fasta.exists():
        fasta.symlink_to(config.GENOME_FASTA)
    pysam.faidx(str(fasta))

    desc = pd.read_csv(config.ANNOTATION_TSV, sep="\t", usecols=["gene_id", "description"])
    desc = dict(zip(desc.gene_id, desc.description.fillna("")))
    names = common_names()

    genes: dict[str, list] = {}
    rows = []
    gid_re = re.compile(r'gene_id "([^"]+)"')
    tid_re = re.compile(r'transcript_id "[^"]+"')
    with open(config.GTF) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gid = gid_re.search(f[8]).group(1)
            f[8] = tid_re.sub(f'transcript_id "{gid}"', f[8]) + f' gene_name "{gid}";'
            rows.append(f)
            if f[2] == "transcript":
                start, end = int(f[3]), int(f[4])
                g = genes.setdefault(gid, [gid, f[0], start, end, f[6], desc.get(gid, ""), names.get(gid, [])])
                g[2], g[3] = min(g[2], start), max(g[3], end)

    rows.sort(key=lambda f: (f[0], int(f[3])))
    gtf = out / "genes.gtf"
    gtf.write_text("".join("\t".join(f) + "\n" for f in rows))
    pysam.tabix_index(str(gtf), preset="gff", force=True)  # -> genes.gtf.gz + .tbi

    (out / "genes.json").write_text(json.dumps(list(genes.values()), separators=(",", ":")))
    named = sum(1 for g in genes.values() if g[6])
    print(f"{len(genes):,} genes ({named:,} with a common name) -> {out}")


if __name__ == "__main__":
    main()
