#!/usr/bin/env python3
"""
Step 1 — collect RefSeq genome accessions for the target taxa.

Reads the species list from config/species.txt and the per-species cap from
config/thresholds.yaml. Writes data/accessions.tsv.

    python scripts/01_collect_accessions.py
    python scripts/01_collect_accessions.py --complete-only   # precision tier
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, ROOT, load_config  # noqa: E402

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon"


def read_species(path: Path) -> list[str]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#")[0].strip()
        if line:
            out.append(line)
    return out


def fetch_species(taxon: str, cap: int, complete_only: bool) -> list[tuple]:
    """Page through the Datasets API until `cap` assemblies are collected."""
    rows: list[tuple] = []
    token = None
    while len(rows) < cap:
        page = min(1000, cap - len(rows))
        url = (f"{API}/{urllib.parse.quote(taxon)}/dataset_report"
               f"?filters.assembly_source=refseq&page_size={page}")
        if complete_only:
            url += "&filters.assembly_level=complete_genome"
        if token:
            url += f"&page_token={token}"
        raw = subprocess.run(["curl", "-s", "--max-time", "120", url],
                             capture_output=True, text=True, check=False).stdout
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError:
            print(f"  ! could not parse response for {taxon}", file=sys.stderr)
            break
        reports = payload.get("reports", [])
        if not reports:
            break
        for rep in reports:
            org = rep["organism"]
            rows.append((
                rep["accession"],
                org["organism_name"],
                org.get("infraspecific_names", {}).get("strain", "-"),
                rep["assembly_info"].get("assembly_level", ""),
            ))
        token = payload.get("next_page_token")
        if not token:
            break
        time.sleep(0.25)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--complete-only", action="store_true",
                    help="restrict to Complete Genome assemblies")
    ap.add_argument("--out", type=Path, default=DATA / "accessions.tsv")
    args = ap.parse_args()

    cfg = load_config()
    cap = cfg["sampling"]["cap_per_species"]
    species = read_species(ROOT / "config" / "species.txt")
    print(f"{len(species)} taxa · cap {cap} per species"
          f"{' · complete genomes only' if args.complete_only else ''}\n")

    seen: set[str] = set()
    rows: list[tuple] = []
    for taxon in species:
        got = [r for r in fetch_species(taxon, cap, args.complete_only)
               if r[0] not in seen]
        seen.update(r[0] for r in got)
        rows.extend(got)
        print(f"  {taxon:<40}{len(got):>5}", flush=True)
        time.sleep(0.25)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join("\t".join(r) for r in rows), encoding="utf-8")

    levels: dict[str, int] = {}
    for r in rows:
        levels[r[3]] = levels.get(r[3], 0) + 1
    print(f"\n{len(rows)} assemblies written to {args.out}")
    for lvl, n in sorted(levels.items(), key=lambda x: -x[1]):
        print(f"  {lvl:<20}{n:>6}")


if __name__ == "__main__":
    main()
