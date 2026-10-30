#!/usr/bin/env python3
"""
Browser-ready reference: genome FASTA + index, annotation GTF, gene search index.
One folder per reference in config.REFS:

  data/ref/<ref>/genome.fasta(.fai)  copy of the FASTA (a copy, not a symlink,
                                     so data/ can move to another disk), indexed
  data/ref/<ref>/genes.gtf.gz(.tbi)  the GTF, sorted, bgzipped, tabix-indexed so
                                     igv.js fetches only the visible region
                                     (the page labels models by gene_id)
  data/ref/<ref>/genes.json          [[gene_id, chr, start, end, strand,
                                       description, [common names]]]

Common names are not a field in either ToxoDB file; they are pulled from the
descriptions: parenthesised symbols in ToxoDB-13 ("SRS29B (SAG1)") and
symbol-like trailing words in ToxoDB-68 ("... AP2XI-2"). ~1,600 genes get one.

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/prep_ref.py [ref ...]
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
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
    with open(config.NAMES_GFF) as fh:
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


def build(name: str, ref: dict, desc: dict[str, str], names: dict[str, list[str]]) -> None:
    out = config.DATA / "ref" / name
    out.mkdir(parents=True, exist_ok=True)

    fasta = out / "genome.fasta"
    if fasta.is_symlink():
        fasta.unlink()
    if not fasta.exists() or fasta.stat().st_size != ref["fasta"].stat().st_size:
        shutil.copyfile(ref["fasta"], fasta)
    pysam.faidx(str(fasta))

    genes: dict[str, list] = {}
    rows = []
    gid_re = re.compile(r'gene_id "([^"]+)"')
    with open(ref["gtf"]) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append(f)
            if f[2] == "transcript":
                gid = gid_re.search(f[8]).group(1)
                start, end = int(f[3]), int(f[4])
                g = genes.setdefault(gid, [gid, f[0], start, end, f[6], desc.get(gid, ""), names.get(gid, [])])
                g[2], g[3] = min(g[2], start), max(g[3], end)

    rows.sort(key=lambda f: (f[0], int(f[3])))
    gtf = out / "genes.gtf"
    gtf.write_text("".join("\t".join(f) + "\n" for f in rows))
    pysam.tabix_index(str(gtf), preset="gff", force=True)  # -> genes.gtf.gz + .tbi

    (out / "genes.json").write_text(json.dumps(list(genes.values()), separators=(",", ":")))
    named = sum(1 for g in genes.values() if g[6])
    print(f"{name}: {len(genes):,} genes ({named:,} with a common name) -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("refs", nargs="*", help=f"default: all of {', '.join(config.REFS)}")
    args = ap.parse_args()
    unknown = set(args.refs) - set(config.REFS)
    if unknown:
        ap.error(f"unknown reference(s): {', '.join(sorted(unknown))}")
    ann = pd.read_csv(config.ANNOTATION_TSV, sep="\t", usecols=["gene_id", "description"])
    desc = dict(zip(ann.gene_id, ann.description.fillna("")))
    names = common_names()
    for name in args.refs or config.REFS:
        build(name, config.REFS[name], desc, names)


if __name__ == "__main__":
    main()
