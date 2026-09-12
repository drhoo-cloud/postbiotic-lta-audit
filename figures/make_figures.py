#!/usr/bin/env python3
"""Regenerate every figure in the manuscript and its supplement.

    python figures/make_figures.py

Figure 1 is composed from two panel scripts; Figures 2, 3, S1, S2 and S3 are
each written by a single script. Output is 600 dpi PNG and LZW-compressed TIFF
at the widths used in the submitted manuscript. figS1 needs
data/tables/TableS9a.csv, which scripts/05_aggregate.py writes.
"""
import os, subprocess, sys
from PIL import Image, ImageDraw, ImageFont
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY = (31, 56, 100)
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"


def run(script):
    r = subprocess.run([sys.executable, os.path.join(HERE, script)],
                       cwd=HERE, capture_output=True, text=True)
    if r.returncode:
        tail = (r.stderr.strip().splitlines() or [""])[-1]
        print(f"  skipped {script}: {tail}")
        return False
    print(f"  {script}")
    return True


def trim(im, pad=12):
    a = np.array(im.convert("RGB"))
    ys = np.where((a.sum(2) < 735).any(1))[0]
    return im.crop((0, max(0, ys.min() - pad), im.width, min(im.height, ys.max() + pad)))


def compose_figure1():
    a = trim(Image.open(os.path.join(HERE, "F1A_chem.png")))
    b = trim(Image.open(os.path.join(HERE, "F1B_strains.png")))
    gap = 66
    out = Image.new("RGB", (a.width, a.height + gap + b.height), "white")
    out.paste(a.convert("RGB"), (0, 0))
    out.paste(b.convert("RGB"), (0, a.height + gap))
    d = ImageDraw.Draw(out)
    f = ImageFont.truetype(BOLD, 46)
    d.text((22, 4), "A", font=f, fill=NAVY)
    d.text((22, a.height + gap - 2), "B", font=f, fill=NAVY)
    out.save(os.path.join(HERE, "Figure1.png"))
    out.save(os.path.join(HERE, "Figure1.tif"), compression="tiff_lzw", dpi=(600, 600))
    print(f"  Figure1: {round(out.width/600*25.4,1)} x {round(out.height/600*25.4,1)} mm")


def flatten():
    """Elsevier wants RGB without an alpha channel; matplotlib leaves one in."""
    import glob
    for f in sorted(glob.glob(os.path.join(HERE, "*.tif"))):
        im = Image.open(f)
        if im.mode == "RGBA":
            bg = Image.new("RGB", im.size, "white")
            bg.paste(im, mask=im.split()[-1])
            bg.save(f, compression="tiff_lzw", dpi=(600, 600))
            print(f"  flattened {os.path.basename(f)}")


if __name__ == "__main__":
    print("panels:")
    run("fig1a_ester_chemistry.py")
    run("fig1b_five_strains.py")
    compose_figure1()
    for s in ("fig2_localization_and_process.py", "fig3_measured_vs_genomic.py",
              "figS1_score_distribution.py", "figS2_contact_residues.py",
              "figS3_mprf.py", "graphical_abstract.py"):
        run(s)
    flatten()
    for f in ("F1A_chem.png", "F1B_strains.png"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            os.remove(p)
