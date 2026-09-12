#!/usr/bin/env python3
"""
Step 3 — layer 2: can the enzyme reach its substrate?

The D-alanine esters of LTA sit on the poly(glycerophosphate) chain outside
the membrane, so an enzyme confined to the cytoplasm cannot act on them.
Each unique DltE protein is submitted to Phobius and assigned to one of three
classes: membrane-anchored, secreted, or cytoplasmic.

    python scripts/03_layer2_localization.py             # all unique proteins
    python scripts/03_layer2_localization.py --limit 200 # one batch

★ RefSeq WP_ accessions are non-redundant, so identical proteins share an
  accession. Screening the unique set therefore covers every genome without
  sampling. Set runtime.contact_email in the config before running: EBI
  requires a real contact address for submitted jobs.

★ This layer is a prediction. Cytoplasmic calls need cell fractionation
  before they can be treated as established.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, fetch_proteins, load_config, require_email  # noqa: E402

POLL_INTERVAL = 2.0
POLL_LIMIT = 40


def unique_proteins() -> dict:
    """Map every distinct DltE accession to the genomes that carry it."""
    src = DATA / "layer1.jsonl"
    if not src.exists():
        sys.exit(f"{src} not found — run 02_layer1_presence.py first")
    index: dict[str, dict] = {}
    with open(src, encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            if "error" in rec or not rec.get("dltE"):
                continue
            acc = rec["dltE"]["id"]
            entry = index.setdefault(acc, {"species": rec["species"], "genomes": 0,
                                           "len": rec["dltE"]["len"],
                                           "score": rec["dltE"]["score"]})
            entry["genomes"] += 1
    return index


def classify(topology: str, tm: int, sp: str) -> str:
    if tm >= 1:
        return "A_TM"
    return "B_SP" if sp == "Y" else "C_CYTO"


def submit(endpoint: str, email: str, fmt: str, acc: str, seq: str) -> dict:
    tag = re.sub(r"[^A-Za-z0-9]", "_", acc)
    fasta = f">{tag}\n{seq}\n"
    job = subprocess.run(
        ["curl", "-s", "--max-time", "60", "-X", "POST", f"{endpoint}/run",
         "--data-urlencode", f"email={email}",
         "--data-urlencode", "stype=protein",
         "--data-urlencode", f"format={fmt}",
         "--data-urlencode", f"sequence={fasta}"],
        capture_output=True, text=True, check=False).stdout.strip()
    if not job.startswith("phobius-"):
        raise RuntimeError(job[:80])

    for _ in range(POLL_LIMIT):
        status = subprocess.run(["curl", "-s", "--max-time", "30",
                                 f"{endpoint}/status/{job}"],
                                capture_output=True, text=True,
                                check=False).stdout.strip()
        if status == "FINISHED":
            break
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise RuntimeError(status)
        time.sleep(POLL_INTERVAL)
    else:
        raise RuntimeError("timed out waiting for Phobius")

    text = subprocess.run(["curl", "-s", "--max-time", "30",
                           f"{endpoint}/result/{job}/out"],
                          capture_output=True, text=True, check=False).stdout
    match = re.search(r"\n\S+\s+(\d+)\s+(\S+)\s+(\S+)", text)
    if not match:
        raise RuntimeError("could not parse Phobius output")
    tm, sp, topology = int(match.group(1)), match.group(2), match.group(3)
    return {"tm": tm, "sp": sp, "topology": topology,
            "cls": classify(topology, tm, sp)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    cfg = load_config()
    email = require_email(cfg)
    L2 = cfg["layer2_localization"]

    index = unique_proteins()
    out_path = DATA / "layer2.jsonl"
    done = {json.loads(l)["id"] for l in open(out_path, encoding="utf-8")
            if json.loads(l).get("cls")} if out_path.exists() else set()
    todo = [a for a in index if a not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(index)} unique proteins covering "
          f"{sum(v['genomes'] for v in index.values())} genomes · "
          f"{len(done)} done · this batch {len(todo)}", flush=True)
    if not todo:
        return

    seqs = fetch_proteins(todo)
    lock, counter, started = threading.Lock(), [0], time.time()
    out = open(out_path, "a", encoding="utf-8")

    def worker(chunk):
        for acc in chunk:
            try:
                rec = {"id": acc, **submit(L2["endpoint"], email,
                                           L2["format"], acc, seqs[acc])}
            except Exception as exc:
                rec = {"id": acc, "error": str(exc)[:80]}
            with lock:
                out.write(json.dumps(rec) + "\n")
                counter[0] += 1
                if counter[0] % 50 == 0:
                    out.flush()
                    el = time.time() - started
                    print(f"  {counter[0]}/{len(todo)}  "
                          f"{el:.0f}s ({el / counter[0]:.2f}s each)", flush=True)

    nthreads = cfg["runtime"]["phobius_threads"]
    threads = [threading.Thread(target=worker, args=(todo[i::nthreads],))
               for i in range(nthreads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    out.close()
    print(f"batch finished · {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
