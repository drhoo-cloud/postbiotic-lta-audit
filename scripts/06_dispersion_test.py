#!/usr/bin/env python3
"""
Step 6 — the dispersion test reported in Section 3 and Table S11d.

Question
    Does removal capability show in the *dispersion* of measured D-alanine
    substitution rather than in its mean? A strain that cannot remove
    D-alanine outside the cell has no enzymatic counterweight to conditions
    that add or strip the ester, so its substitution level should be the
    less constrained of the two.

Inclusion (all four must hold; applied before any value was inspected)
    1. a measured D-alanine substitution value in Table S1
    2. a deposited genome assembly
    3. butanol extraction — the largest methodological confounder held fixed
    4. membership of one of the thirty species carried through the
       localization screen (Table S9d)

Range rule
    Table S1 reports some values as ranges. The midpoint is used, applied to
    every qualifying entry and not selectively.

Test
    Two-sided permutation test on absolute deviations from the group median,
    20,000 permutations of the pooled deviations. Dispersion, not location,
    is the quantity of interest; the group means are reported alongside only
    to show that they do not differ.

    python scripts/06_dispersion_test.py
    python scripts/06_dispersion_test.py --include-outside-screen

Result at the time of writing: n = 4 with predicted removal (SD 15.7,
CV 28%) against n = 8 without (SD 24.2, CV 62%); variance ratio 2.37,
two-sided p = 0.70. The comparison is limited by the size of the measured
record, not by the size of any effect.
See the note at the foot of this file on the two Apilactobacillus species
that fall outside the screen.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import statistics as st
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
S1 = ROOT / "data" / "tables" / "TableS1_structure_map.csv"
FIG2 = ROOT / "figures" / "fig2_data.json"

N_PERM = 20_000
SEED = 0

# Genus abbreviations used in fig2_data.json, expanded to match Table S1.
ABBR = {
    "Lac.": "Lactococcus", "Li.": "Limosilactobacillus", "Lb.": "Lactobacillus",
    "Lg.": "Ligilactobacillus", "Lt.": "Latilactobacillus", "P.": "Pediococcus",
    "Lc.": "Lacticaseibacillus", "Lv.": "Levilactobacillus",
    "Lp.": "Lactiplantibacillus", "Leu.": "Leuconostoc", "Ap.": "Apilactobacillus",
}
# Table S1 abbreviates the genus on repeat mentions; resolve against the
# first full mention of the same species earlier in the table.
CONTINUED = {
    "L. plantarum": "Lactiplantibacillus plantarum",
    "L. rhamnosus": "Lacticaseibacillus rhamnosus",
    "L. delbrueckii": "Lactobacillus delbrueckii",
}
NON_LAB = {"Staphylococcus", "Listeria", "Streptococcus"}

# Species with a measured value that were not carried through the screen
# (fewer than 25 genomes available), listed so the exclusion is explicit.
OUTSIDE_SCREEN = {"Apilactobacillus kosoi", "Apilactobacillus apinorum"}


def removal_class() -> dict[str, str]:
    """Species -> 'removal' | 'no removal', from the localization screen."""
    data = json.loads(FIG2.read_text())
    out = {}
    for short, v in data.items():
        abbr, epithet = short.split(" ", 1)
        reaches = (v["counts"][0] + v["counts"][1]) / v["n"]
        out[f"{ABBR[abbr]} {epithet}"] = "removal" if reaches >= 0.5 else "no removal"
    return out


def midpoint(cell: str):
    """'40-55' -> 47.5 ; '42' -> 42.0 ; 'NR' -> None."""
    cell = (cell or "").strip()
    if not cell or cell.upper() in {"NR", "ND", "NA", "-"}:
        return None
    nums = [float(t) for t in re.findall(r"\d+\.?\d*", cell)]
    return sum(nums) / len(nums) if nums else None


def load(include_outside: bool):
    cls = removal_class()
    rows = []
    for r in csv.DictReader(S1.open(encoding="utf-8-sig")):
        strain = (r.get("Strain") or "").strip()
        if not strain or strain.startswith("Note"):
            continue
        species = " ".join(strain.split()[:2])
        species = CONTINUED.get(species, species)
        if species.split()[0] in NON_LAB:
            continue
        value = midpoint(r.get("D-alanine (%)"))
        if value is None:
            continue
        if "BuOH" not in (r.get("Extraction method") or ""):
            continue
        if not (r.get("Genome accession") or "").strip() or \
                "no assembly" in (r.get("Genome accession") or ""):
            continue
        group = cls.get(species)
        if group is None:
            if include_outside and species in OUTSIDE_SCREEN:
                group = "no removal"
            else:
                continue
        rows.append({"strain": strain, "species": species,
                     "group": group, "value": value,
                     "source": (r.get("D-alanine (%)") or "").strip()})
    return rows


def permutation_p(a: np.ndarray, b: np.ndarray, n_perm: int, seed: int) -> float:
    """Two-sided test on the difference in mean absolute deviation."""
    rng = np.random.default_rng(seed)
    da = np.abs(a - np.median(a))
    db = np.abs(b - np.median(b))
    obs = abs(db.mean() - da.mean())
    pool = np.concatenate([da, db])
    n1 = len(da)
    hits = 0
    for _ in range(n_perm):
        p = rng.permutation(pool)
        if abs(p[n1:].mean() - p[:n1].mean()) >= obs:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--include-outside-screen", action="store_true",
                    help="add the two Apilactobacillus species that were not "
                         "carried through the screen (reported only as a "
                         "sensitivity check; see the note in this file)")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "tables" / "TableS11d.csv")
    args = ap.parse_args()

    if not S1.exists():
        sys.exit(f"{S1} not found")
    rows = load(args.include_outside_screen)
    if not rows:
        sys.exit("no qualifying strains — check Table S1")

    rows.sort(key=lambda r: (r["group"] != "removal", -r["value"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["strain", "species", "predicted_class",
                    "value_used_pct", "source_range"])
        for r in rows:
            w.writerow([r["strain"], r["species"], r["group"],
                        f"{r['value']:.1f}", r["source"]])

    print(f"{len(rows)} strains from "
          f"{len({r['species'] for r in rows})} species\n")
    print(f"{'strain':<34}{'class':<12}{'value':>7}   source")
    for r in rows:
        print(f"{r['strain'][:33]:<34}{r['group']:<12}{r['value']:>7.1f}   {r['source']}")

    groups = {}
    for name in ("removal", "no removal"):
        v = np.array([r["value"] for r in rows if r["group"] == name], float)
        groups[name] = v
        print(f"\n{name}: n={len(v)} mean={v.mean():.1f} "
              f"SD={v.std(ddof=1):.1f} CV={v.std(ddof=1) / v.mean() * 100:.0f}% "
              f"range {v.min():.0f}-{v.max():.0f}")

    a, b = groups["removal"], groups["no removal"]
    ratio = b.var(ddof=1) / a.var(ddof=1)
    p = permutation_p(a, b, N_PERM, SEED)
    print(f"\nvariance ratio (no removal / removal) = {ratio:.2f}")
    print(f"permutation p (two-sided, {N_PERM:,} permutations, seed {SEED}) = {p:.3f}")
    print(f"\nwritten to {args.out}")


# ---------------------------------------------------------------------------
# Note on the two excluded Apilactobacillus species
#
# Ap. kosoi 10HT and Ap. apinorum JCM 30765T both carry a deposited genome and
# a measured substitution of zero. They are excluded because fewer than 25
# genomes were available for either species, so neither was carried through
# the localization screen. Including them raises the variance ratio to 2.96
# and leaves the test null (p = 0.12).
#
# They would not belong on this axis in any case. Neither species lacks the
# removal enzyme alone: both
# lack the dlt operon altogether, and their lipoteichoic acid carries L-lysine
# in place of D-alanine. The substituent whose dispersion is being compared is
# not present in them, so they do not sit on this axis.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
