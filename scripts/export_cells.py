#!/usr/bin/env python3
"""
Per-cell table for the page, one columnar JSON per line.

2D UMAP and UMI totals come from the line's combined.h5ad (UT and IAA share
one embedding). Clusters, phase, pseudotime and stage calls come from the
Sc-x2 coords table. Leiden clusters are per condition, so the page labels them
"<condition> · <n>".

Output: data/<line>/cells.json
  {line, conditions: [...], bams: {condition: "<file>.bam"},
   cb, cond (index into conditions), x, y, cluster, phase, pseudotime,
   brady_ees, stage10, umis}

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/export_cells.py AP2XI-2_AP2XII-1
"""
from __future__ import annotations

import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd

import config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("line", choices=config.LINES)
    args = ap.parse_args()
    line = config.LINES[args.line]

    a = ad.read_h5ad(line["h5ad"], backed="r")
    umap = pd.DataFrame(a.obsm["X_umap"][:, :2], index=a.obs_names, columns=["x", "y"])
    umap["umis"] = a.obs["total_counts"].to_numpy()

    coords = pd.read_csv(config.COORDS)
    coords = coords[coords.dataset == line["coords_dataset"]].set_index("barcode")
    df = coords.join(umap, how="inner")
    missing = len(coords) - len(df)
    if missing:
        print(f"warning: {missing} cells in the coords table have no 2D UMAP")

    conditions = list(line["samples"])
    df = df[df.condition.isin(conditions)]
    out = {
        "line": line["label"],
        "conditions": conditions,
        "bams": {c: f"{config.sample_id(c)}.bam" for c in conditions},
        "cb": [config.bam_cb(b) for b in df.index],
        "cond": df.condition.map(conditions.index).tolist(),
        "x": np.round(df.x.astype(float), 3).tolist(),
        "y": np.round(df.y.astype(float), 3).tolist(),
        "cluster": df.leiden_cluster.astype(int).tolist(),
        "phase": df.phase_call.fillna("").tolist(),
        "pseudotime": np.round(df.pseudotime.astype(float), 4).fillna(-1).tolist(),
        "brady_ees": df.brady_ees_call.fillna("").tolist(),
        "stage10": df.ees_10stage_call.fillna("").tolist(),
        "umis": df.umis.astype(int).tolist(),
    }
    dest = config.DATA / args.line / "cells.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, separators=(",", ":")))
    print(f"{len(df):,} cells -> {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
