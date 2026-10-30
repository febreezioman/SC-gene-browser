#!/usr/bin/env python3
"""
Shrink a Cell Ranger BAM to browser size and tag reads with their cluster.

For each condition of a line:
  1. samtools keeps mapped reads whose CB is a called cell in the coords table
     and drops UMI duplicates (0x400), so each molecule shows once. Antisense,
     intronic and multimapped reads are kept.
  2. pysam adds CL:Z:<leiden cluster> so igv.js can group reads by cluster.
  3. Sort is preserved; the result is indexed.

The page filters by the CB tag itself, so no read-group rewrite is needed.

Output: data/<line>/<condition>.bam (+ .bai)

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/prep_bam.py AP2XI-2_AP2XII-1
       [--keep-duplicates] [--threads 8]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time

import pandas as pd
import pysam

import config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("line", choices=config.LINES)
    ap.add_argument("--condition", help="only this condition (default: all)")
    ap.add_argument("--keep-duplicates", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    line = config.LINES[args.line]
    cells = pd.read_csv(config.COORDS, usecols=["dataset", "condition", "barcode", "leiden_cluster"])
    cells = cells[cells.dataset == line["coords_dataset"]]

    outdir = config.DATA / args.line
    outdir.mkdir(parents=True, exist_ok=True)

    for condition, run_dir in line["samples"].items():
        if args.condition and condition != args.condition:
            continue
        sub = cells[cells.condition == condition]
        cluster = dict(zip(sub.barcode.map(config.bam_cb), sub.leiden_cluster.astype(str)))
        sid = config.sample_id(condition)
        bc_file = outdir / f"{sid}.barcodes.txt"
        bc_file.write_text("\n".join(cluster) + "\n")

        src = run_dir / "outs" / "possorted_genome_bam.bam"
        tmp = outdir / f"{sid}.tmp.bam"
        out = outdir / f"{sid}.bam"
        flags = "0x4" if args.keep_duplicates else "0x404"
        t0 = time.time()
        print(f"[{condition}] {len(cluster):,} cells; filtering {src}", flush=True)
        subprocess.run(
            ["samtools", "view", "-b", "-@", str(args.threads), "-F", flags,
             "-D", f"CB:{bc_file}", "-o", str(tmp), str(src)],
            check=True,
        )

        n = 0
        with pysam.AlignmentFile(tmp, "rb", threads=args.threads) as fin, \
             pysam.AlignmentFile(out, "wb", template=fin, threads=args.threads) as fout:
            for read in fin:
                read.set_tag("CL", cluster[read.get_tag("CB")], "Z")
                fout.write(read)
                n += 1
        tmp.unlink()
        pysam.index(str(out), "-@", str(args.threads))
        size = out.stat().st_size / 1e9
        print(f"[{condition}] {n:,} reads -> {out} ({size:.2f} GB, {time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
