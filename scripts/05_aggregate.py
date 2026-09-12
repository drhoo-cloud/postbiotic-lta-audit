#!/usr/bin/env python3
"""
Step 5 — aggregate the three layers into the supplementary tables.

Writes per-genome calls, per-species prevalence with Wilson intervals, and
the localization breakdown. Controls are tabulated separately from the lactic
acid bacterial totals, as specified in the configuration.

    python scripts/05_aggregate.py
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, load_config, normalise_species, wilson  # noqa: E402


def load_layers() -> tuple[list, dict]:
    l1 = DATA / "layer1.jsonl"
    if not l1.exists():
        sys.exit(f"{l1} not found — run 02_layer1_presence.py first")
    genomes = []
    with open(l1, encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            if "error" not in rec:
                rec["sp"] = normalise_species(rec["species"])
                genomes.append(rec)

    calls: dict[str, str] = {}
    l2 = DATA / "layer2.jsonl"
    if l2.exists():
        with open(l2, encoding="utf-8") as fh:
            for line in fh:
                rec = json.loads(line)
                if rec.get("cls"):
                    calls[rec["id"]] = rec["cls"]
    else:
        print("! layer2.jsonl absent — localization columns will be empty",
              file=sys.stderr)
    return genomes, calls


def classify(rec: dict, calls: dict) -> str:
    if not rec.get("dltA"):
        return "E_no_dlt"
    if not rec.get("dltE"):
        return "D_dltE_absent"
    return calls.get(rec["dltE"]["id"], "unassigned")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=DATA / "tables")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cfg = load_config()
    controls = set(cfg["sampling"]["controls"])
    scores = sorted(set(cfg["layer1_presence"]["sensitivity_bitscores"]))
    primary = cfg["layer1_presence"]["ortholog_call"]["min_bitscore"]

    genomes, calls = load_layers()
    for rec in genomes:
        rec["loc"] = classify(rec, calls)

    # ── per-genome table
    cols = ["accession", "species", "strain", "assembly_level", "n_protein",
            "dltA", "dltE", "dltE_len", "dltE_score",
            *[f"dltE_s{s}" for s in scores],
            "mprF_long", "mprF_short", "localization"]
    with open(args.outdir / "S9a_per_genome.csv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in genomes:
            e = r.get("dltE")
            w.writerow([r["acc"], r["species"], r.get("strain", ""), r["level"],
                        r["nprot"], int(bool(r.get("dltA"))), int(bool(e)),
                        e["len"] if e else "", e["score"] if e else "",
                        *[int(bool(r.get("dltE" if s == primary else f"dltE_s{s}")))
                          for s in scores],
                        int(bool(r.get("mprF_long"))),
                        int(bool(r.get("mprF_short"))), r["loc"]])

    # ── per-species table
    by_sp = collections.defaultdict(list)
    for r in genomes:
        by_sp[r["sp"]].append(r)
    with open(args.outdir / "S9b_per_species.csv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "is_control", "n_total", "n_complete", "dltA",
                    "dltE", "prevalence", "wilson_lo", "wilson_hi",
                    *[f"dltE_s{s}" for s in scores],
                    "A_TM", "B_SP", "C_CYTO", "D_absent", "E_no_dlt",
                    "mprF_long_pct", "mprF_short_pct"])
        for sp in sorted(by_sp):
            rs = by_sp[sp]
            n = len(rs)
            e = sum(1 for r in rs if r.get("dltE"))
            lo, hi = wilson(e, n)
            loc = collections.Counter(r["loc"] for r in rs)
            w.writerow([sp, int(sp in controls), n,
                        sum(1 for r in rs if r["level"] == "Complete Genome"),
                        sum(1 for r in rs if r.get("dltA")), e,
                        f"{e / n:.4f}", f"{lo:.4f}", f"{hi:.4f}",
                        *[sum(1 for r in rs
                              if r.get("dltE" if s == primary else f"dltE_s{s}"))
                          for s in scores],
                        loc["A_TM"], loc["B_SP"], loc["C_CYTO"],
                        loc["D_dltE_absent"], loc["E_no_dlt"],
                        f"{sum(1 for r in rs if r.get('mprF_long')) / n:.3f}",
                        f"{sum(1 for r in rs if r.get('mprF_short')) / n:.3f}"])

    # ── headline numbers
    lab = [r for r in genomes if r["sp"] not in controls]
    n = len(lab)
    print(f"{len(genomes)} genomes · {len(by_sp)} species "
          f"({n} lactic acid bacteria, {len(genomes) - n} controls)\n")
    counts = collections.Counter(r["loc"] for r in lab)
    labels = {"A_TM": "membrane-anchored", "B_SP": "secreted",
              "C_CYTO": "cytoplasmic (predicted)", "D_dltE_absent": "DltE absent",
              "E_no_dlt": "dlt operon absent", "unassigned": "unassigned"}
    for key, name in labels.items():
        if counts[key]:
            lo, hi = wilson(counts[key], n)
            print(f"  {name:<26}{counts[key]:>6}  {counts[key] / n:>6.1%}"
                  f"  ({lo:.1%}-{hi:.1%})")
    no_access = counts["C_CYTO"] + counts["D_dltE_absent"] + counts["E_no_dlt"]
    lo, hi = wilson(no_access, n)
    print(f"\n  no extracytoplasmic removal capacity: "
          f"{no_access}/{n} = {no_access / n:.1%} ({lo:.1%}-{hi:.1%})")
    # ── figure input: species with enough genomes to plot
    fig_rows = []
    for sp, rs in by_sp.items():
        n_sp = len(rs)
        if n_sp < 25:
            continue
        loc = collections.Counter(r["loc"] for r in rs)
        e = sum(1 for r in rs if r.get("dltE"))
        lo, hi = wilson(e, n_sp)
        fig_rows.append({
            "sp": sp, "n": n_sp, "ctrl": int(sp in controls),
            "A": loc["A_TM"], "B": loc["B_SP"], "C": loc["C_CYTO"],
            "Dd": loc["D_dltE_absent"], "E": loc["E_no_dlt"],
            "prev": e / n_sp, "lo": lo, "hi": hi,
            "p90": sum(1 for r in rs if r.get(f"dltE_s{scores[0]}")) / n_sp,
            "p130": sum(1 for r in rs if r.get(f"dltE_s{scores[-1]}")) / n_sp,
        })
    fig_rows.sort(key=lambda r: (-(r["C"] + r["Dd"] + r["E"]) / r["n"], -r["n"]))
    (DATA / "figdata.json").write_text(json.dumps(fig_rows))
    print(f"  figure input: {len(fig_rows)} species -> {DATA / 'figdata.json'}")

    print(f"\ntables written to {args.outdir}")


if __name__ == "__main__":
    main()
