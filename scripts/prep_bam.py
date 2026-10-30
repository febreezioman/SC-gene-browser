#!/usr/bin/env python3
"""
Shrink a line's BAMs to browser size and tag reads with their group.

Reads come from the cells in data/<line>/cells.json (run export_cells.py
first), and each kept read gets CL:Z:<value of the line's cl_field> so igv.js
can group reads by it. Antisense, intronic, intergenic and spliced reads are
all kept; only reads from non-cells and duplicate reads of a molecule go.

  cellranger  samtools keeps the cells' reads and drops reads flagged 0x400
              (Cell Ranger marks all but one read per molecule); pysam adds CL.
  starsolo    STARsolo flags no duplicates, so one read is kept per molecule:
              per (CB, UB, GX gene) when the read is assigned to a gene, else
              per (CB, UB, strand, start). Reads without a UMI or with a
              secondary/supplementary flag are dropped. Contigs run in
              parallel into config.SCRATCH and are concatenated in header
              order (stays sorted).

Output: data/<line>/<condition>.bam (+ .bai)

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/prep_bam.py <line>
       [--condition C] [--threads 12]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import pysam

import config

_CL: dict[str, str] = {}  # CB -> CL value, shared with forked workers


def cellranger(src: Path, out: Path, cl: dict[str, str], threads: int) -> int:
    bc_file = out.with_suffix(".barcodes.txt")
    bc_file.write_text("\n".join(cl) + "\n")
    tmp = out.with_suffix(".tmp.bam")
    subprocess.run(["samtools", "view", "-b", "-@", str(threads), "-F", "0x404",
                    "-D", f"CB:{bc_file}", "-o", str(tmp), str(src)], check=True)
    n = 0
    with pysam.AlignmentFile(tmp, "rb", threads=threads) as fin, \
         pysam.AlignmentFile(out, "wb", template=fin, threads=threads) as fout:
        for read in fin:
            read.set_tag("CL", cl[read.get_tag("CB")], "Z")
            fout.write(read)
            n += 1
    tmp.unlink()
    return n


def _dedup_contig(job: tuple[str, str, str]) -> tuple[str, int, int, int]:
    contig, src, part = job
    kept = dups = no_umi = 0
    seen = set()
    with pysam.AlignmentFile(src) as fin, pysam.AlignmentFile(part, "wb", template=fin) as fout:
        for r in fin.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            cl = _CL.get(r.get_tag("CB"))
            if cl is None:
                continue
            ub = r.get_tag("UB")
            if ub == "-":
                no_umi += 1
                continue
            gx = r.get_tag("GX") if r.has_tag("GX") else "-"
            key = (r.get_tag("CB"), ub, gx) if gx != "-" else (r.get_tag("CB"), ub, r.is_reverse, r.reference_start)
            if key in seen:
                dups += 1
                continue
            seen.add(key)
            r.set_tag("CL", cl, "Z")
            fout.write(r)
            kept += 1
    return contig, kept, dups, no_umi


def starsolo(src: Path, out: Path, cl: dict[str, str], threads: int) -> int:
    _CL.clear()
    _CL.update(cl)
    parts = config.SCRATCH / f"{out.parent.name}_{out.stem}.parts"
    parts.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(src) as f:
        stats = {s.contig: s.mapped for s in f.get_index_statistics()}
        contigs = [c for c in f.references if stats.get(c, 0) > 0]
    jobs = [(c, str(src), str(parts / f"{i:05d}.bam")) for i, c in enumerate(contigs)]
    totals = [0, 0, 0]
    # Biggest contigs first so the pool is not left waiting on one at the end.
    with Pool(threads) as pool:
        for contig, kept, dups, no_umi in pool.imap_unordered(_dedup_contig, sorted(jobs, key=lambda j: -stats[j[0]])):
            totals = [totals[0] + kept, totals[1] + dups, totals[2] + no_umi]
            if stats[contig] > 5_000_000:
                print(f"  {contig}: kept {kept:,}, duplicates {dups:,}", flush=True)
    print(f"  kept {totals[0]:,} reads; dropped {totals[1]:,} duplicates, {totals[2]:,} without UMI", flush=True)
    subprocess.run(["samtools", "cat", "-o", str(out)] + [j[2] for j in jobs], check=True)
    for j in jobs:
        Path(j[2]).unlink()
    parts.rmdir()
    return totals[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("line", choices=config.LINES)
    ap.add_argument("--condition", help="only this condition (default: all)")
    ap.add_argument("--threads", type=int, default=12)
    args = ap.parse_args()

    line = config.LINES[args.line]
    outdir = config.DATA / args.line
    cells = json.loads((outdir / "cells.json").read_text())
    values = cells["values"][cells["cl_field"]]

    for k, (condition, src) in enumerate(line["samples"].items()):
        if args.condition and condition != args.condition:
            continue
        cl = {cb: values[i] for i, cb in enumerate(cells["cb"]) if cells["cond"][i] == k}
        out = outdir / f"{config.sample_id(condition)}.bam"
        t0 = time.time()
        print(f"[{condition}] {len(cl):,} cells; {src}", flush=True)
        n = (cellranger if line["kind"] == "cellranger" else starsolo)(Path(src), out, cl, args.threads)
        pysam.index(str(out), "-@", str(args.threads))
        print(f"[{condition}] {n:,} reads -> {out} ({out.stat().st_size / 1e9:.2f} GB, {time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
