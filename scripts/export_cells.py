#!/usr/bin/env python3
"""
Per-cell table for the page, one columnar JSON per line, plus data/lines.json
(the line switcher's list).

Each line declares its own fields, so the page builds its colour menu, groups
and hover text from the data:
  fields: [{key, label, type: "cat"|"num", order?, prefix?}]   values: {key: [...]}
`cl_field` names the field that prep_bam.py writes into each read's CL tag.

  AP2 lines (Cell Ranger): 2D UMAP and UMI totals from the line's combined.h5ad
    (conditions share one embedding); clusters, phase, pseudotime and stage
    calls from the Sc-x2 coords table. Leiden clusters are per condition.
  AlRubaye 2026 cat gut (STARsolo): UMAP, Leiden 0.4 and the cell types
    transferred from the published atlas, from Sc-x2 scripts/59 and 62.

Output: data/<line>/cells.json, data/lines.json

Usage: ~/projects/Sc-x2/envs/sc/bin/python scripts/export_cells.py <line>
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import config

PHASES = ["G1a", "G1b", "S", "M", "C"]
BRADY_EES = ["Tachyzoite", "Bradyzoite", "Bradyzoite/EES1-3", "EES1-3"]
STAGES10 = ["Tachyzoites", "BRP1+", "Merozoite A", "Merozoite B", "Merozoite C",
            "Trans. Merozoites", "Pre-gamete", "Female", "Male-1", "Male-2"]
# Order used in the AlRubaye figures (Sc-x2 scripts/62).
ALRUBAYE_TYPES = ["Merozoite-A", "Merozoite-B", "Merozoite-C", "Transitioning-Merozoites",
                  "BRP1+", "pre-gametes", "Female", "Male-1", "Tachyzoites"]


def cellranger_line(line: dict) -> tuple[pd.DataFrame, list[dict], str]:
    import anndata as ad

    a = ad.read_h5ad(line["h5ad"], backed="r")
    umap = pd.DataFrame(a.obsm["X_umap"][:, :2], index=a.obs_names, columns=["x", "y"])
    umap["umis"] = a.obs["total_counts"].to_numpy()
    coords = pd.read_csv(config.COORDS)
    coords = coords[coords.dataset == line["coords_dataset"]].set_index("barcode")
    df = coords.join(umap, how="inner")
    if len(df) < len(coords):
        print(f"warning: {len(coords) - len(df)} cells in the coords table have no 2D UMAP")
    df = pd.DataFrame({
        "cb": [config.bam_cb(b) for b in df.index],
        "condition": df.condition, "x": df.x, "y": df.y,
        "cluster": df.leiden_cluster.astype(int).astype(str),
        "phase": df.phase_call, "stage10": df.ees_10stage_call, "brady_ees": df.brady_ees_call,
        "pseudotime": df.pseudotime, "umis": df.umis,
    }).reset_index(drop=True)
    fields = [
        {"key": "cluster", "label": "Leiden cluster", "type": "cat", "prefix": "Cluster "},
        {"key": "phase", "label": "Cell-cycle phase", "type": "cat", "order": PHASES},
        {"key": "stage10", "label": "Stage (10-stage call)", "type": "cat", "order": STAGES10},
        {"key": "brady_ees", "label": "Bradyzoite / EES call", "type": "cat", "order": BRADY_EES},
        {"key": "pseudotime", "label": "Pseudotime", "type": "num"},
        {"key": "umis", "label": "UMIs per cell", "type": "num"},
    ]
    return df, fields, "cluster"


def alrubaye_line(line: dict) -> tuple[pd.DataFrame, list[dict], str]:
    cells = pd.read_csv(config.ALRUBAYE_CELLS / "cells.tsv", sep="\t", dtype={"guide_target": str})
    umap = pd.read_csv(config.ALRUBAYE_CELLS / "viz" / "umap_coords.tsv.gz", sep="\t", usecols=["barcode", "umap1", "umap2", "leiden"])
    df = cells.merge(umap, on="barcode", how="inner")
    (condition,) = line["samples"]
    guides = df.guide_target.fillna("")
    df = pd.DataFrame({
        "cb": df.barcode, "condition": condition, "x": df.umap1, "y": df.umap2,
        "celltype": df.label, "leiden": df.leiden.astype(str),
        "guide": guides.where(guides != "", "none"),
        "r": df.r, "umis": df.umis,
    })
    guide_order = ["none"] + [g for g in df.guide.value_counts().index if g != "none"]
    fields = [
        {"key": "celltype", "label": "Cell type (transferred from atlas)", "type": "cat", "order": ALRUBAYE_TYPES},
        {"key": "leiden", "label": "Leiden cluster (0.4)", "type": "cat", "prefix": "Cluster "},
        {"key": "guide", "label": "Guide target (single call)", "type": "cat", "order": guide_order},
        {"key": "r", "label": "Label confidence (r)", "type": "num"},
        {"key": "umis", "label": "Parasite UMIs per cell", "type": "num"},
    ]
    return df, fields, "celltype"


def column(values: pd.Series, kind: str) -> list:
    if kind == "num":
        v = pd.to_numeric(values, errors="coerce").astype(float).round(4)
        return [None if np.isnan(x) else x for x in v]
    return values.fillna("").astype(str).tolist()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("line", choices=config.LINES)
    args = ap.parse_args()
    line = config.LINES[args.line]

    df, fields, cl_field = (cellranger_line if line["kind"] == "cellranger" else alrubaye_line)(line)
    conditions = list(line["samples"])
    df = df[df.condition.isin(conditions)]

    out = {
        "line": line["label"],
        "ref": line["ref"],
        "conditions": conditions,
        "bams": {c: f"{config.sample_id(c)}.bam" for c in conditions},
        "cl_field": cl_field,
        "cb": df.cb.tolist(),
        "cond": df.condition.map(conditions.index).tolist(),
        "x": np.round(df.x.astype(float), 3).tolist(),
        "y": np.round(df.y.astype(float), 3).tolist(),
        "fields": fields,
        "values": {f["key"]: column(df[f["key"]], f["type"]) for f in fields},
    }
    dest = config.DATA / args.line / "cells.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, separators=(",", ":")))
    print(f"{len(df):,} cells -> {dest} ({dest.stat().st_size / 1e6:.1f} MB)")

    index_file = config.DATA / "lines.json"
    lines = json.loads(index_file.read_text()) if index_file.exists() else {}
    lines[args.line] = {"label": line["label"], "ref": line["ref"], "cells": len(df), "conditions": conditions}
    index_file.write_text(json.dumps(lines, indent=1))


if __name__ == "__main__":
    main()
