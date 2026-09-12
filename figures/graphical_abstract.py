"""Graphical abstract - Elsevier specification.

Light tinted panels, one accent per panel, no heavy greys. Reads left to
right: what a postbiotic must carry, why sequence cannot supply the
missing one, what does set it, and the four entries that close it. Each
panel ends on a concept strip rather than a number.

Typography follows the Elsevier artwork note: Arial-metric sans only, and
sized to stay legible when ScienceDirect scales the image into a
500 x 200 px window.

190 x 76 mm at 600 dpi = 4488 x 1795 px, ratio 2.50 against the required
minimum of 1328 x 531 px.

    python figures/graphical_abstract.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyBboxPatch, Rectangle, Circle
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
F = 'Liberation Sans'          # Arial-metric; rasterised, so nothing is embedded
plt.rcParams.update({'font.family': F, 'font.sans-serif': [F],
                     'mathtext.fontset': 'custom', 'mathtext.rm': F,
                     'mathtext.it': F + ':italic', 'mathtext.default': 'regular'})

INK = '#243B5E'          # headings and numbers
BLACK = '#111111'        # concept strips
BODY = '#3C4657'         # running text
SUB = '#7C8798'          # sub-labels
LINE = '#D9E1EC'
STRIP = '#E7EFF9'        # concept strip
ARROW = '#8FB4DE'
AMBER_F = '#FDF3E2'; AMBER_E = '#E0B463'; AMBER_T = '#9A5B1E'

TINT = [('#E6F0FA', '#F8FBFE', '#C7DAF0'),    # panel 1  blue
        ('#EFEBF8', '#FBFAFE', '#D5CCEC'),    # panel 2  lavender
        ('#E6F3EC', '#F8FCFA', '#C6E2D2'),    # panel 3  mint
        ('#FBEBEC', '#FEF9F9', '#EFC9CC')]    # panel 4  rose

MM = 1 / 25.4
fig = plt.figure(figsize=(190 * MM, 76 * MM), dpi=600)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis('off')

# ---------------- banner ----------------
ax.add_patch(Rectangle((0.0, 33.8), 100.0, 6.2, fc='#DCE8F6', ec='none', zorder=1))
ax.text(50.0, 36.7,
        'STANDARDIZING POSTBIOTICS   \u00b7   WHAT TO SPECIFY FOR LIPOTEICHOIC ACID   '
        '\u00b7   WHAT SEQUENCE CANNOT SUPPLY',
        ha='center', va='center', fontsize=9.0, color=INK,
        fontweight='bold', zorder=4)

TOP, BOT = 31.6, 1.4
HEAD_H = 4.4
CONC_H = 4.6


def panel(x, w, i, title):
    head, body, edge = TINT[i]
    ax.add_patch(FancyBboxPatch((x, BOT), w, TOP - BOT,
                                boxstyle='round,pad=0,rounding_size=1.2',
                                fc=body, ec=edge, lw=1.2, zorder=2))
    ax.add_patch(Rectangle((x + 0.9, TOP - HEAD_H - 0.9), w - 1.8, HEAD_H,
                           fc=head, ec='none', zorder=3))
    ax.text(x + w / 2, TOP - 0.9 - HEAD_H / 2, title, ha='center', va='center',
            fontsize=8.2, color=INK, fontweight='bold', zorder=4)


def concept(x, w, text):
    ax.add_patch(Rectangle((x + 0.8, BOT + 0.9), w - 1.6, CONC_H,
                           fc=STRIP, ec='none', zorder=3))
    ax.text(x + w / 2, BOT + 0.9 + CONC_H / 2, text, ha='center', va='center',
            fontsize=7.0, color=BLACK, fontweight='bold', linespacing=1.45, zorder=4)


def arrow(x, dx=2.6):
    ax.add_patch(FancyArrow(x, (TOP + BOT) / 2, dx, 0, width=0.8, head_width=2.4,
                            head_length=1.4, length_includes_head=True,
                            fc=ARROW, ec='none', zorder=3))


def stat(xc, y, big, small):
    ax.text(xc, y, big, ha='center', va='center', fontsize=8.7, color=INK,
            fontweight='bold', zorder=4)
    ax.text(xc, y - 2.7, small, ha='center', va='center', fontsize=7.2,
            color=SUB, zorder=4)


# ---------------- 1 ----------------
X1, W1 = 1.0, 23.5
panel(X1, W1, 0, 'A POSTBIOTIC MUST CARRY')
ax.text(X1 + 4.6, TOP - 10.2, 'which organism, and how many', ha='left', va='center',
        fontsize=7.0, color=BODY, zorder=4)
ax.plot([X1 + 2.1, X1 + 2.8, X1 + 3.9], [TOP - 10.2, TOP - 11.1, TOP - 9.1],
        color=ARROW, lw=1.9, solid_capstyle='round', zorder=4)
ax.add_patch(FancyBboxPatch((X1 + 1.6, TOP - 19.4), W1 - 3.2, 5.6,
                            boxstyle='round,pad=0,rounding_size=0.7',
                            fc=AMBER_F, ec=AMBER_E, lw=1.1, zorder=3))
ax.text(X1 + W1 / 2, TOP - 16.6, 'how much of the\nactive component',
        ha='center', va='center', fontsize=7.6, color=AMBER_T,
        fontweight='bold', linespacing=1.4, zorder=4)
ax.text(X1 + 4.6, TOP - 22.4, 'how it was inactivated', ha='left', va='center',
        fontsize=7.0, color=BODY, zorder=4)
ax.plot([X1 + 2.1, X1 + 2.8, X1 + 3.9], [TOP - 22.4, TOP - 23.3, TOP - 21.3],
        color=ARROW, lw=1.9, solid_capstyle='round', zorder=4)
concept(X1, W1, 'ONE OF THE THREE HAS\nNO WORKED INSTANCE')

arrow(25.1)

# ---------------- 2 ----------------
X2, W2 = 28.7, 21.0
panel(X2, W2, 1, 'SEQUENCE CANNOT')
gy = TOP - 7.6
ax.add_patch(FancyArrow(X2 + 2.4, gy, 8.6, 0, width=1.9, head_width=1.9,
                        head_length=1.1, length_includes_head=True,
                        fc='#C9D7EA', ec='none', zorder=3))
ax.text(X2 + 6.5, gy, '$dlt$ A\u2013D', ha='center', va='center', fontsize=7.0,
        color=INK, fontweight='bold', zorder=4)
ax.add_patch(FancyArrow(X2 + 12.0, gy, 6.2, 0, width=1.9, head_width=1.9,
                        head_length=1.1, length_includes_head=True,
                        fc='#E3E3E8', ec='none', zorder=3))
ax.text(X2 + 14.9, gy, '$dltE$', ha='center', va='center', fontsize=7.0,
        color=SUB, fontweight='bold', zorder=4)
ax.text(X2 + 6.5, gy - 2.6, 'in every panel', ha='center', va='center',
        fontsize=6.8, color=BODY, zorder=4)
ax.text(X2 + 14.9, gy - 2.6, 'in none', ha='center', va='center',
        fontsize=6.8, color=AMBER_T, fontweight='bold', zorder=4)
stat(X2 + W2 / 2, TOP - 13.6, u'3.2\u00d7', 'apart, with identical proteins')
ax.plot([X2 + 3.2, X2 + W2 - 3.2], [TOP - 18.2] * 2, color=LINE, lw=1.0, zorder=3)
stat(X2 + W2 / 2, TOP - 20.1, '56.6%', 'of 7,404 genomes, no removal')
concept(X2, W2, 'THE GENOME GIVES THE\nMEANS, NOT THE VALUE')

arrow(50.3)

# ---------------- 3 ----------------
X3, W3 = 53.9, 21.0
panel(X3, W3, 2, 'THE PROCESS DOES')
fy = TOP - 7.6
STEPS = ['ferment', 'inactivate', 'extract']
sw, gap = 5.4, 1.1
sx = X3 + (W3 - (3 * sw + 2 * gap)) / 2
for i, t in enumerate(STEPS):
    x = sx + i * (sw + gap)
    hot = (i == 1)
    ax.add_patch(FancyBboxPatch((x, fy - 1.5), sw, 3.0,
                                boxstyle='round,pad=0,rounding_size=0.5',
                                fc='#CFE3D6' if hot else '#EDF5F0',
                                ec='none', zorder=3))
    ax.text(x + sw / 2, fy, t, ha='center', va='center', fontsize=6.4,
            color=INK if hot else BODY, fontweight='bold' if hot else 'normal',
            zorder=4)
    if i < 2:
        ax.plot([x + sw + 0.12, x + sw + gap - 0.12], [fy, fy], color=LINE,
                lw=1.1, zorder=3)
ax.text(X3 + W3 / 2, fy - 3.2, 'every step moves the value', ha='center',
        va='center', fontsize=6.8, color=SUB, zorder=4)
stat(X3 + W3 / 2, TOP - 13.6, u'10.7\u00d7', u'culture pH 6 \u2192 8, one strain')
ax.plot([X3 + 3.2, X3 + W3 - 3.2], [TOP - 18.2] * 2, color=LINE, lw=1.0, zorder=3)
stat(X3 + W3 / 2, TOP - 20.1, u'60\u201369%', 'of the ester, stripped by heat')
concept(X3, W3, 'THE PROCESS CARRIES\nTHE VALUE')

arrow(75.5)

# ---------------- 4 ----------------
X4, W4 = 79.1, 19.9
panel(X4, W4, 3, 'SO STATE IT')
FOUR = ['substitution level', 'chain length', 'lipoprotein control', 'process record']
for i, t in enumerate(FOUR):
    y = TOP - 9.6 - i * 4.4
    ax.add_patch(Circle((X4 + 3.0, y), 1.05, fc='#B8555C', ec='none', zorder=3))
    ax.text(X4 + 3.0, y, str(i + 1), ha='center', va='center', fontsize=6.4,
            color='white', fontweight='bold', zorder=4)
    ax.text(X4 + 5.0, y, t, ha='left', va='center', fontsize=7.6, color=BODY, zorder=4)
concept(X4, W4, 'MEASURED, NOT INFERRED')

fig.savefig(os.path.join(HERE, 'GA.png'), dpi=600, facecolor='white')
p = os.path.join(HERE, 'GA.tif')
fig.savefig(p, dpi=600, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
im = Image.open(p)
if im.mode == 'RGBA':
    bg = Image.new('RGB', im.size, 'white')
    bg.paste(im, mask=im.split()[-1])
    bg.save(p, compression='tiff_lzw', dpi=(600, 600))
print('graphical abstract written:', Image.open(p).size, Image.open(p).mode)
