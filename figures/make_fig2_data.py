#!/usr/bin/env python3
"""Regenerate figures/fig2_data.json from data/tables/TableS9d.csv.

The JSON is a convenience copy of the per-species localization table so the
figure scripts run without a CSV parser. It is derived, never edited by hand.

    python figures/make_fig2_data.py
"""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data", "tables", "TableS9d.csv")
DST = os.path.join(HERE, "fig2_data.json")

ABBR = {"Lactococcus": "Lac.", "Limosilactobacillus": "Li.", "Lactobacillus": "Lb.",
        "Ligilactobacillus": "Lg.", "Latilactobacillus": "Lt.", "Pediococcus": "P.",
        "Lacticaseibacillus": "Lc.", "Levilactobacillus": "Lv.",
        "Lactiplantibacillus": "Lp.", "Leuconostoc": "Leu.", "Apilactobacillus": "Ap."}
COLS = ["A membrane-anchored", "B secreted", "C cytoplasmic",
        "D DltE absent", "E dlt operon absent"]

out = {}
with open(SRC, encoding="utf-8-sig") as fh:
    for r in csv.DictReader(fh):
        name = (r.get("Species") or "").strip()
        if not name or not (r.get("Genomes") or "").strip().isdigit():
            continue
        genus, epithet = name.split()[0], " ".join(name.split()[1:])
        if genus not in ABBR:
            continue
        counts = [int(r[c]) for c in COLS]
        n = int(r["Genomes"])
        if sum(counts) != n:
            raise SystemExit(f"{name}: counts sum to {sum(counts)}, not {n}")
        out[f"{ABBR[genus]} {epithet}"] = {
            "n": n, "counts": counts,
            "frac": [round(c / n * 100, 2) for c in counts]}

with open(DST, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"{len(out)} species written to {DST}")
