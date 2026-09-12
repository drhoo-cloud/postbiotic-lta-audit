"""
Shared helpers: configuration, paths, reference sequences, statistics.

All thresholds live in config/thresholds.yaml. Nothing in this repository
should hard-code a cut-off; if you find one, it is a bug.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "thresholds.yaml"
DATA = Path(os.environ.get("LTA_DATA_DIR", ROOT / "data"))
DATA.mkdir(parents=True, exist_ok=True)
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# ── configuration ─────────────────────────────────────────────────────
def load_config(path: Path | str = CONFIG) -> dict:
    with open(path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    """Fail loudly rather than silently screening with the wrong criteria."""
    call = cfg["layer1_presence"]["ortholog_call"]
    lo, hi = call["length_range"]
    if not 0 < lo < hi:
        raise ValueError(f"invalid length_range: {call['length_range']}")
    if call["min_bitscore"] not in cfg["layer1_presence"]["sensitivity_bitscores"]:
        raise ValueError(
            "min_bitscore must appear in sensitivity_bitscores so that the "
            "primary call is reported alongside its sensitivity analysis"
        )
    if cfg["statistics"]["pass_threshold"] is not None:
        raise ValueError(
            "pass_threshold must stay null. Reporting a distribution rather "
            "than a pass/fail call is part of the design, not an oversight."
        )


def require_email(cfg: dict) -> str:
    """EBI asks for a real contact address on every submitted job."""
    email = cfg.get("runtime", {}).get("contact_email")
    if not email or "@" not in str(email):
        sys.exit(
            "Set runtime.contact_email in config/thresholds.yaml to your own\n"
            "address before running the Phobius step. EBI requires a real\n"
            "contact for submitted jobs; please do not use a placeholder."
        )
    return email


# ── FASTA ─────────────────────────────────────────────────────────────
def read_fasta(path: Path | str) -> tuple[dict, dict]:
    """Return (accession -> sequence, accession -> product annotation)."""
    seqs, ann, cur = {}, {}, None
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                acc, _, rest = line[1:].partition(" ")
                cur = acc
                seqs[cur] = ""
                ann[cur] = rest.split(" [")[0].strip()
            elif cur:
                seqs[cur] += line.strip()
    return seqs, ann


def fetch_proteins(accessions: list[str], batch: int = 150) -> dict:
    """Fetch protein sequences from NCBI by accession.

    Reference sequences are retrieved this way rather than extracted from a
    downloaded genome, so the pipeline runs without a pre-staged assembly.
    """
    import time

    out: dict[str, str] = {}
    for i in range(0, len(accessions), batch):
        chunk = ",".join(accessions[i : i + batch])
        text = subprocess.run(
            ["curl", "-s", "-G", f"{NCBI_EUTILS}/efetch.fcgi", "--data",
             f"db=protein&rettype=fasta&retmode=text&id={chunk}"],
            capture_output=True, text=True, check=False,
        ).stdout
        cur = None
        for line in text.splitlines():
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = ""
            elif cur:
                out[cur] += line.strip()
        time.sleep(0.4)
    missing = [a for a in accessions if a not in out]
    if missing:
        raise RuntimeError(f"could not fetch reference sequences: {missing}")
    return out


def reference_sequences(cfg: dict) -> dict:
    """Reference proteins named in the config, cached under data/."""
    cache = DATA / "reference_proteins.json"
    wanted = {k: v["accession"] for k, v in cfg["reference"].items()
              if isinstance(v, dict) and "accession" in v}
    if cache.exists():
        cached = json.loads(cache.read_text())
        if set(cached) == set(wanted.values()):
            return {name: cached[acc] for name, acc in wanted.items()}
    seqs = fetch_proteins(sorted(set(wanted.values())))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(seqs))
    return {name: seqs[acc] for name, acc in wanted.items()}


# ── statistics ────────────────────────────────────────────────────────
def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval. Returns (0, 0) for an empty sample."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def normalise_species(name: str) -> str:
    for token in ("subsp.", "="):
        name = name.split(token)[0]
    return " ".join(name.split()[:2])


# ── resumable JSONL ───────────────────────────────────────────────────
def completed(path: Path | str, key: str = "acc") -> set:
    """Accessions already processed successfully, for resuming a run."""
    done: set = set()
    if not Path(path).exists():
        return done
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "error" not in rec and key in rec:
                done.add(rec[key])
    return done
