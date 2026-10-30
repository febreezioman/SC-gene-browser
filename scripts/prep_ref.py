#!/usr/bin/env python3
"""
Browser-ready reference: genome FASTA + index, annotation GTF, gene search index.

  data/ref/genome.fasta(.fai)  symlink to the Sc-x2 FASTA, indexed here
  data/ref/genes.gtf.gz(.tbi)  the Cell Ranger GTF with gene_name = gene_id
                               (igv.js labels models by gene ID), sorted,
                               bgzipped and tabix-indexed so igv.js fetches
                               only the visible region
  data/ref/genes.json          [[gene_id, chr, start, end, strand, description]]
                               for the page's gene search (ToxoDB-68 descriptions)

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/prep_ref.py
"""
from __future__ import annotations

import json
import re

import pandas as pd
import pysam

import config


def main() -> None:
    out = config.DATA / "ref"
    out.mkdir(parents=True, exist_ok=True)

    fasta = out / "genome.fasta"
    if not fasta.exists():
        fasta.symlink_to(config.GENOME_FASTA)
    pysam.faidx(str(fasta))

    desc = pd.read_csv(config.ANNOTATION_TSV, sep="\t", usecols=["gene_id", "description"])
    desc = dict(zip(desc.gene_id, desc.description.fillna("")))

    genes: dict[str, list] = {}
    rows = []
    gid_re = re.compile(r'gene_id "([^"]+)"')
    with open(config.GTF) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gid = gid_re.search(f[8]).group(1)
            f[8] += f' gene_name "{gid}";'
            rows.append(f)
            if f[2] == "transcript":
                start, end = int(f[3]), int(f[4])
                g = genes.setdefault(gid, [gid, f[0], start, end, f[6], desc.get(gid, "")])
                g[2], g[3] = min(g[2], start), max(g[3], end)

    rows.sort(key=lambda f: (f[0], int(f[3])))
    gtf = out / "genes.gtf"
    gtf.write_text("".join("\t".join(f) + "\n" for f in rows))
    pysam.tabix_index(str(gtf), preset="gff", force=True)  # -> genes.gtf.gz + .tbi

    (out / "genes.json").write_text(json.dumps(list(genes.values()), separators=(",", ":")))
    print(f"{len(genes):,} genes -> {out}")


if __name__ == "__main__":
    main()
