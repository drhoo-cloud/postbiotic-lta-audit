#!/usr/bin/env python3
"""
lipo_detect.py -- lipobox-based lipoprotein prediction for Gram-positive proteomes.

Predictor used for the a1 analysis (SignalP/LipoP unavailable offline).

Call rule (all four must hold):
  1. n-region charge : >=1 K/R within residues 2..8
  2. Cys position    : invariant Cys at residue 15..40 (1-based)
  3. lipobox         : [LVIFMAGWST][ASGTVILMF][GAS] immediately before that Cys
  4. h-region        : mean Kyte-Doolittle >= 1.5 over the 8 residues ending at
                       lipobox position -4, and no D/E/K/R/P in the lipobox itself
"""
import re, sys

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,'X':0.0,'U':0.0,'B':-3.5,'Z':-3.5}

LIPOBOX = re.compile(r'[LVIFMAGWST][ASGTVILMF][GAS]C')
CYS_MIN, CYS_MAX = 15, 40
WINDOW = 60


def kd_mean(seq):
    return sum(KD.get(c, 0.0) for c in seq) / len(seq) if seq else 0.0


def predict(seq):
    """Return (is_lipoprotein, detail_dict)."""
    s = seq.upper().replace('*', '')
    head = s[:WINDOW]
    d = {'cys_pos': None, 'lipobox': None, 'h_kd': None, 'n_charge': None,
         'reason': None}

    n_reg = head[1:8]
    d['n_charge'] = n_reg.count('K') + n_reg.count('R')
    if d['n_charge'] < 1:
        d['reason'] = 'no positive charge in n-region'
        return False, d

    for m in LIPOBOX.finditer(head):
        cys = m.start() + 4               # 1-based position of the Cys
        if not (CYS_MIN <= cys <= CYS_MAX):
            continue
        box = m.group(0)[:3]
        if any(c in box for c in 'DEKRP'):
            continue
        h_start = max(0, m.start() - 8)
        h_seq = head[h_start:m.start()]
        h = kd_mean(h_seq)
        if h < 1.5:
            d.update(cys_pos=cys, lipobox=box, h_kd=round(h, 2),
                     reason='h-region not hydrophobic enough')
            continue
        d.update(cys_pos=cys, lipobox=box, h_kd=round(h, 2), reason='lipobox')
        return True, d

    if d['reason'] is None:
        d['reason'] = 'no lipobox with Cys in window'
    return False, d


def read_fasta(path):
    name, buf = None, []
    fh = sys.stdin if path == '-' else open(path)
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield name, ''.join(buf)
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name:
        yield name, ''.join(buf)
    if fh is not sys.stdin:
        fh.close()


if __name__ == '__main__':
    n_tot = n_lipo = 0
    for name, seq in read_fasta(sys.argv[1]):
        ok, d = predict(seq)
        n_tot += 1
        n_lipo += ok
        if '-v' in sys.argv or ok:
            print(f"{'LIPO' if ok else '----'}\t{name.split()[0]}\t"
                  f"cys={d['cys_pos']}\tbox={d['lipobox']}\tkd={d['h_kd']}\t{d['reason']}")
    print(f"# total={n_tot} lipoproteins={n_lipo} "
          f"density={n_lipo/n_tot*100:.2f}%" if n_tot else "# empty", file=sys.stderr)
