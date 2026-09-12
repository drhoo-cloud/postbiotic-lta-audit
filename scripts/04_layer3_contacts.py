#!/usr/bin/env python3
"""
Step 4 — layer 3: how tightly does the enzyme hold the substrate?

PDB 8AKH is the DltE extracellular domain soaked with LTA, with
sn-glycerol-3-phosphate — the LTA backbone repeat unit — bound in the active
site. Residues within the configured cut-off of that ligand are mapped across
every unique DltE protein.

    python scripts/04_layer3_contacts.py                # uses cached contacts
    python scripts/04_layer3_contacts.py --recompute    # re-derive from PDB

The split between catalytic/scaffold and substrate-contact residues follows
the published motif assignment and was not re-grouped after seeing results.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import Background, Builder

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, fetch_proteins, load_config, reference_sequences, wilson  # noqa: E402

ABC = Alphabet.amino()
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
       "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
       "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PDBE = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}.cif"


def _name(x) -> str:
    return x.decode() if isinstance(x, bytes) else x


def parse_residue(label: str) -> tuple[str, int]:
    """'Ser128' -> ('S', 128)."""
    return AA3[label[:3]], int(label[3:])


def derive_contacts(cfg: dict) -> list[tuple[str, float]]:
    """Recompute ligand contacts from the deposited structure."""
    import gemmi

    struct_cfg = cfg["reference"]["structure"]
    pdb = struct_cfg["pdb"].lower()
    cif = DATA / f"{pdb}.cif"
    if not cif.exists():
        cif.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sL", "--max-time", "120", "-o", str(cif),
                        PDBE.format(pdb=pdb)], check=False)

    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_alternative_conformations()
    model = st[0]
    ligand = struct_cfg["ligand"]
    skip = {"HOH", "GOL", ligand}
    cutoff = cfg["layer3_contacts"]["contact_cutoff_angstrom"]

    ns = gemmi.NeighborSearch(st, 5.0).populate()
    nearest: dict[tuple[str, int], float] = {}
    for chain in model:
        for res in chain:
            if res.name != ligand:
                continue
            for atom in res:
                for mark in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                    cra = mark.to_cra(model)
                    if cra.residue.name in skip:
                        continue
                    dist = cra.atom.pos.dist(atom.pos)
                    key = (cra.residue.name.capitalize(), cra.residue.seqid.num)
                    if key not in nearest or dist < nearest[key]:
                        nearest[key] = dist
    return sorted(((f"{n}{i}", round(d, 1)) for (n, i), d in nearest.items()),
                  key=lambda x: x[1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--recompute", action="store_true",
                    help="re-derive contact residues from the PDB entry")
    args = ap.parse_args()

    cfg = load_config()
    L3 = cfg["layer3_contacts"]
    groups = {r: "catalytic_scaffold" for r in L3["catalytic_scaffold"]}
    groups.update({r: "substrate_contact" for r in L3["substrate_contact"]})

    if args.recompute:
        found = derive_contacts(cfg)
        print(f"contacts within {L3['contact_cutoff_angstrom']} A of the ligand:")
        for label, dist in found:
            mark = "" if label in groups else "   ← not in config"
            print(f"  {label:<10}{dist:>5} A{mark}")
        missing = [r for r in groups if r not in dict(found)]
        if missing:
            print(f"\n! configured residues absent from the structure: {missing}")
        print()

    # ── unique DltE proteins from layer 1
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
    print(f"{len(index)} unique proteins · "
          f"{sum(v['genomes'] for v in index.values())} genomes")

    seqs = fetch_proteins(sorted(index))
    ref = reference_sequences(cfg)["dltE"]

    # ── align everything to a single-sequence HMM built from the reference
    bg = Background(ABC)
    hmm, _, _ = Builder(ABC).build(
        TextSequence(name=b"ref", sequence=ref).digitize(ABC), bg)
    names = sorted(seqs) + ["__REF__"]
    block = DigitalSequenceBlock(ABC, [
        TextSequence(name=n.encode(),
                     sequence=(ref if n == "__REF__" else seqs[n])).digitize(ABC)
        for n in names])
    msa = pyhmmer.hmmer.hmmalign(hmm, block)
    rows = {_name(n): (a if isinstance(a, str) else a.decode())
            for n, a in zip(msa.names, msa.alignment)}

    # map reference numbering onto alignment columns
    column, pos = {}, 0
    for idx, ch in enumerate(rows["__REF__"]):
        if ch not in "-.":
            pos += 1
            column[pos] = idx

    labels = L3["catalytic_scaffold"] + L3["substrate_contact"]
    counts = {lab: collections.Counter() for lab in labels}
    weighted = {lab: collections.Counter() for lab in labels}
    out_rows = []
    for acc in sorted(seqs):
        aln = rows[acc]
        meta = index[acc]
        row = [acc, meta["species"], meta["len"], meta["score"], meta["genomes"]]
        for lab in labels:
            _, num = parse_residue(lab)
            col = column.get(num)
            value = (aln[col].upper()
                     if col is not None and col < len(aln) and aln[col] not in "-."
                     else "-")
            counts[lab][value] += 1
            weighted[lab][value] += meta["genomes"]
            row.append(value)
        out_rows.append(row)

    dest = DATA / "layer3_contacts.csv"
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein_id", "species", "length", "hmm_score",
                         "n_genomes"] + labels)
        writer.writerows(out_rows)

    n_prot = len(out_rows)
    n_gen = sum(r[4] for r in out_rows)
    print(f"\n{'group':<20}{'residue':<10}{'proteins':>10}{'genomes':>10}"
          f"{'Wilson 95% CI':>20}")
    print("-" * 72)
    for lab in labels:
        wt, _ = parse_residue(lab)
        lo, hi = wilson(weighted[lab][wt], n_gen)
        print(f"{groups[lab]:<20}{lab:<10}"
              f"{counts[lab][wt] / n_prot:>9.1%}{weighted[lab][wt] / n_gen:>10.1%}"
              f"{f'{lo:.1%}-{hi:.1%}':>20}")
    for group in ("catalytic_scaffold", "substrate_contact"):
        members = [l for l in labels if groups[l] == group]
        mean = sum(counts[l][parse_residue(l)[0]] for l in members) / (len(members) * n_prot)
        print(f"{group + ' mean':<31}{mean:>9.1%}")
    print(f"\nwritten to {dest}")


if __name__ == "__main__":
    main()
