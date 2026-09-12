#!/usr/bin/env python3
"""
Step 2 — layer 1: is the D-alanine removal enzyme present?

Downloads each proteome, searches it with HMMs built from the reference
proteins, and records the call at every sensitivity threshold. Results are
appended to data/layer1.jsonl, so an interrupted run resumes where it stopped.

    python scripts/02_layer1_presence.py                 # whole list
    python scripts/02_layer1_presence.py --limit 500     # one batch

Every cut-off comes from config/thresholds.yaml. None is hard-coded here.
"""
from __future__ import annotations

import argparse
import json
import os
import queue
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile, TextSequence
from pyhmmer.plan7 import Background, Builder

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, completed, load_config, read_fasta, reference_sequences  # noqa: E402

ABC = Alphabet.amino()
DL = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
      "{acc}/download?include_annotation_type=PROT_FASTA")


def _name(x) -> str:
    return x.decode() if isinstance(x, bytes) else x


def build_hmms(cfg: dict):
    refs = reference_sequences(cfg)
    bg = Background(ABC)
    builder = Builder(ABC)
    hmms = []
    for label, seq in refs.items():
        hmm, _, _ = builder.build(
            TextSequence(name=label.encode(), sequence=seq).digitize(ABC), bg)
        hmms.append(hmm)
    return hmms


def downloader(items, out_q, tmp: Path, retries: int = 4):
    for acc, species, strain, level in items:
        zpath, fpath, ok = tmp / f"{acc}.zip", tmp / f"{acc}.faa", False
        for attempt in range(retries):
            try:
                subprocess.run(
                    ["curl", "-sL", "--max-time", "90", "--retry", "2",
                     "-o", str(zpath), DL.format(acc=acc)], check=False)
                with zipfile.ZipFile(zpath) as zf:
                    for member in zf.namelist():
                        if member.endswith("protein.faa"):
                            fpath.write_bytes(zf.read(member))
                            ok = True
                            break
                zpath.unlink(missing_ok=True)
                if ok:
                    break
            except Exception:
                zpath.unlink(missing_ok=True)
            time.sleep(1.5 * (attempt + 1))
        out_q.put((acc, species, strain, level, fpath if ok else None))


def screen(cfg: dict, limit: int | None) -> None:
    acc_file = DATA / "accessions.tsv"
    if not acc_file.exists():
        sys.exit(f"{acc_file} not found — run 01_collect_accessions.py first")
    rows = [ln.split("\t") for ln in acc_file.read_text().splitlines() if ln.strip()]

    out_path = DATA / "layer1.jsonl"
    done = completed(out_path)
    todo = [r for r in rows if r[0] not in done]
    if limit:
        todo = todo[:limit]
    print(f"{len(done)}/{len(rows)} already screened · this batch {len(todo)}",
          flush=True)
    if not todo:
        return

    L1 = cfg["layer1_presence"]
    call = L1["ortholog_call"]
    lo_len, hi_len = call["length_range"]
    evalue = float(L1["evalue"])
    scores = sorted(set(L1["sensitivity_bitscores"] + [call["min_bitscore"]]))
    hmms = build_hmms(cfg)

    tmp = Path(os.environ.get("TMPDIR", "/tmp"))
    q: queue.Queue = queue.Queue(maxsize=24)
    nthreads = cfg["runtime"]["download_threads"]
    for i in range(nthreads):
        threading.Thread(target=downloader,
                         args=(todo[i::nthreads], q, tmp), daemon=True).start()

    started, ok = time.time(), 0
    with open(out_path, "a", encoding="utf-8") as out:
        for i in range(len(todo)):
            acc, species, strain, level, fpath = q.get()
            try:
                if fpath is None:
                    raise RuntimeError("download failed")
                _, ann = read_fasta(fpath)
                with SequenceFile(fpath, digital=True, alphabet=ABC) as fh:
                    targets = list(fh)
                lengths = {_name(s.name): len(s.sequence) for s in targets}
                block = DigitalSequenceBlock(ABC, targets)

                hits: dict[str, list] = {}
                for res in pyhmmer.hmmer.hmmsearch(
                        hmms, block, E=evalue, cpus=cfg["runtime"]["hmmer_cpus"]):
                    key = _name(res.query.name)
                    found = [(_name(h.name), lengths[_name(h.name)], round(h.score, 1))
                             for h in res if h.evalue < evalue]
                    hits[key] = sorted(found, key=lambda x: -x[2])[:5]

                def pick(label, lo, hi, min_score):
                    for name, length, score in hits.get(label, []):
                        if lo <= length <= hi and score >= min_score:
                            return {"id": name, "len": length,
                                    "score": score, "ann": ann.get(name, "")}
                    return None

                rec = {"acc": acc, "species": species, "strain": strain,
                       "level": level, "nprot": len(targets),
                       "dltA": pick("dltA", 400, 600, 150),
                       "mprF_long": pick("mprF_long", 700, 1050, 300),
                       "mprF_short": pick("mprF_short", 260, 430, 120)}
                for s in scores:
                    tag = "dltE" if s == call["min_bitscore"] else f"dltE_s{s}"
                    rec[tag] = pick("dltE", lo_len, hi_len, s)
                out.write(json.dumps(rec, ensure_ascii=False) + "\n")
                ok += 1
            except Exception as exc:
                out.write(json.dumps({"acc": acc, "species": species,
                                      "level": level,
                                      "error": str(exc)[:80]}) + "\n")
            finally:
                if fpath is not None:
                    Path(fpath).unlink(missing_ok=True)
            if (i + 1) % 200 == 0:
                out.flush()
                el = time.time() - started
                print(f"  {i + 1}/{len(todo)}  ok={ok}  "
                      f"{el:.0f}s ({el / (i + 1):.2f}s per genome)", flush=True)
    print(f"batch finished · {ok} screened · {time.time() - started:.0f}s")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after N genomes (the run is resumable)")
    args = ap.parse_args()
    screen(load_config(), args.limit)


if __name__ == "__main__":
    main()
