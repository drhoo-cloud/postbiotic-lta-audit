"""Figure 2 — what sets the value.

Panel A: predicted DltE localization across the thirty screened species,
from data/tables/TableS9d.csv via fig2_data.json.

Panel B: how far each variable moves the D-alanine substitution level,
expressed as fold reduction from the least-stressed condition in the same
study, so that series reported in different units can be set side by side.
Source values are printed against each row.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np, json

HERE = os.path.dirname(os.path.abspath(__file__))
F = 'Liberation Sans'
plt.rcParams.update({'font.family': F, 'font.sans-serif': [F],
                     'mathtext.fontset': 'custom', 'mathtext.rm': F,
                     'mathtext.it': F + ':italic', 'mathtext.bf': F + ':bold',
                     'mathtext.default': 'regular',
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})

NAVY = '#1F3864'; BLUE = '#5B8FC9'; RED = '#C00000'; PINK = '#E8A0A0'
GREY = '#BFBFBF'; TXT = '#333333'; AX = '#9A9A9A'
COLS = [NAVY, BLUE, RED, PINK, GREY]
KEYS = ['Membrane-anchored', 'Secreted', 'Cytoplasmic (predicted)',
        'DltE absent', '$\\it{dlt}$ operon absent']
MM = 1 / 25.4

D = json.load(open(os.path.join(HERE, 'fig2_data.json')))
# same ordering as scripts/05_aggregate.py: no-access fraction first, then size
ORDER = sorted(D, key=lambda k: (-(sum(D[k]['counts'][2:]) / D[k]['n']), -D[k]['n']))
TOTAL = [22.7, 20.8, 21.8, 33.3, 1.4]

fig = plt.figure(figsize=(174 * MM, 168 * MM), dpi=600)

# ================= PANEL A : localization =================
L, R = 0.300, 0.925
ax = fig.add_axes([L, 0.470, R - L, 0.460])
ax.set_xlim(0, 100)
N = len(ORDER)
ax.set_ylim(N + 1.9, -1.0)
ax.axis('off')


def stack(y, frac, h=0.62, lw=.35):
    x = 0.0
    for f, c in zip(frac, COLS):
        if f <= 0:
            x += f; continue
        ax.add_patch(Rectangle((x, y - h / 2), f, h, fc=c, ec='white', lw=lw, zorder=3))
        x += f


YB = -0.35
stack(YB, TOTAL, h=1.15, lw=.5)
ax.text(-1.6, YB, 'All lactic acid bacteria', ha='right', va='center',
        fontsize=7.0, color=NAVY, fontweight='bold')
ax.text(101.6, YB, '7,404', ha='left', va='center', fontsize=7.0,
        color=NAVY, fontweight='bold')
xs = 0.0
for v, c in zip(TOTAL, COLS):
    if v > 5:
        ax.text(xs + v / 2, YB, f'{v}%', ha='center', va='center', fontsize=6.6,
                color='white', fontweight='bold', zorder=5)
    xs += v
bx0 = TOTAL[0] + TOTAL[1]
ax.plot([bx0, bx0, 100, 100], [-1.20, -1.42, -1.42, -1.20], color=RED, lw=.8,
        zorder=4, clip_on=False)
ax.text((bx0 + 100) / 2, -1.75, 'no extracytoplasmic access   56.6% (55.4\u201357.7%)',
        ha='center', va='center', fontsize=6.6, color=RED, fontweight='bold',
        clip_on=False)
ax.plot([bx0] * 2, [YB - 0.60, YB + 0.60], color='#2B2B2B', lw=1.2, zorder=6)
ax.plot([-2, 102], [0.42, 0.42], color='#D8D8D8', lw=.6, clip_on=False)

full = [i for i, nm in enumerate(ORDER) if sum(D[nm]['counts'][2:]) == D[nm]['n']]
Y0, Y1 = min(full) + 1.6, max(full) + 1.6
ax.add_patch(Rectangle((-30, Y0 - 0.52), 141, Y1 - Y0 + 1.04, fc='#F6EBEB',
                       ec='none', zorder=0, clip_on=False))
ax.text(108.5, (Y0 + Y1) / 2, 'no removal in any genome', ha='center', va='center',
        fontsize=6.4, color=RED, fontweight='bold', rotation=90, clip_on=False)

ax.text(101.6, 0.95, '$\\it{n}$', ha='left', va='center', fontsize=6.4, color=TXT)
for i, nm in enumerate(ORDER):
    y = i + 1.6
    stack(y, D[nm]['frac'])
    bnd = D[nm]['frac'][0] + D[nm]['frac'][1]
    if 0.4 < bnd < 99.6:
        ax.plot([bnd, bnd], [y - 0.36, y + 0.36], color='#2B2B2B', lw=1.0, zorder=6)
    ax.text(-1.6, y, nm, ha='right', va='center', fontsize=5.4, color=TXT, style='italic')
    ax.text(101.6, y, str(D[nm]['n']), ha='left', va='center', fontsize=6.0, color=TXT)

for xv in (0, 25, 50, 75, 100):
    ax.plot([xv, xv], [N + 1.15, N + 1.52], color='#4D4D4D', lw=1.0, clip_on=False)
    ax.text(xv, N + 2.05, f'{xv}%' if xv == 100 else str(xv), ha='center', va='center',
            fontsize=6.2, color=TXT, clip_on=False)
ax.plot([0, 100], [N + 1.15, N + 1.15], color='#4D4D4D', lw=1.2, clip_on=False)
ax.plot([0, 0], [1.05, N + 1.15], color='#4D4D4D', lw=1.2, clip_on=False)
fig.text((L + R) / 2, 0.452, 'Genomes per species (%)', ha='center', va='center',
         fontsize=6.8, color=TXT)

handles = [Rectangle((0, 0), 1, 1, fc=c, ec='none') for c in COLS]
leg = ax.legend(handles, KEYS, loc='lower center', bbox_to_anchor=(0.5, 1.045), ncol=5,
                frameon=False, fontsize=6.2, handlelength=1.0, handleheight=0.85,
                handletextpad=0.42, columnspacing=1.2, borderpad=0)
for t in leg.get_texts():
    t.set_color(TXT)

# ================= PANEL B : what moves the value =================
# (label, from, to, unit note, fold, source, group)
ROWS = [
 ('Growth pH 6.07 \u2192 8.10', '0.75 \u2192 0.07 mol/mol P', 10.7, 'MacArthur & Archibald 1984', 'process'),
 ('Extraction butanol \u2192 phenol', '69 \u2192 ~25% of units', 2.8, 'R\u00e4is\u00e4nen et al. 2007', 'method'),
 ('NaCl 2 \u2192 100 g/L', '0.73 \u2192 0.30 mol/mol Gro', 2.4, 'Fischer & R\u00f6sel 1980', 'process'),
 ('Heat 52 \u00b0C, 15 min', '0.67 \u2192 0.30 nmol/mg', 2.2, 'Hurst et al. 1975', 'process'),
 ('Five strains, one method', '55 \u2192 17% of units', 3.2, 'Shiraishi et al. 2025', 'strain'),
 ('Two strains, same medium', '0.73 vs 0.62 mol/mol Gro', 1.18, 'Fischer & R\u00f6sel 1980', 'strain'),
]
CG = {'process': NAVY, 'method': RED, 'strain': GREY}
bx = fig.add_axes([0.300, 0.075, 0.470, 0.300])
y = np.arange(len(ROWS))
bx.barh(y, [r[2] for r in ROWS], height=0.58,
        color=[CG[r[4]] for r in ROWS], zorder=3)
for i, r in enumerate(ROWS):
    bx.text(r[2] + 0.22, i, f'{r[2]:.1f}\u00d7' if r[2] >= 2 else f'{r[2]:.2f}\u00d7',
            ha='left', va='center', fontsize=6.4, color=CG[r[4]], fontweight='bold')
bx.axvline(1.0, color='#4D4D4D', lw=.8, ls=(0, (3, 2)), zorder=2)
bx.set_yticks(y)
bx.set_yticklabels([r[0] for r in ROWS], fontsize=6.2, color=TXT)
bx.set_ylim(len(ROWS) - 0.4, -0.7)
bx.set_xlim(0, 12.6)
bx.set_xticks([0, 2, 4, 6, 8, 10, 12])
bx.set_xlabel('Fold reduction in D-alanine substitution, within one study',
              fontsize=6.8, color=TXT, labelpad=3)
bx.tick_params(labelsize=6.2, colors=TXT, length=3.0, width=1.0, color='#4D4D4D', pad=2)
for sp in ('top', 'right'):
    bx.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    bx.spines[sp].set_color('#4D4D4D'); bx.spines[sp].set_linewidth(1.0)
bx.grid(axis='x', color='#EDEDED', lw=.6, zorder=0)
bx.set_axisbelow(True)

for i, r in enumerate(ROWS):
    fig.text(0.790, 0.075 + 0.300 * (1 - (i + 0.5) / len(ROWS)),
             f'{r[1]}\n{r[3]}', ha='left', va='center', fontsize=5.4,
             color='#6E6E6E', linespacing=1.35)

lg = [Rectangle((0, 0), 1, 1, fc=CG[k], ec='none') for k in ('process', 'method', 'strain')]
lb = ['Fermentation or inactivation', 'Extraction method', 'Strain identity']
lg2 = bx.legend(lg, lb, loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3,
                frameon=False, fontsize=6.2, handlelength=1.0, handleheight=0.85,
                handletextpad=0.42, columnspacing=1.2, borderpad=0)
for t in lg2.get_texts():
    t.set_color(TXT)

fig.text(0.020, 0.985, 'A', ha='left', va='top', fontsize=9, fontweight='bold', color='black')
fig.text(0.020, 0.420, 'B', ha='left', va='top', fontsize=9, fontweight='bold', color='black')

fig.savefig(os.path.join(HERE, 'Figure2.png'), dpi=600, facecolor='white')
fig.savefig(os.path.join(HERE, 'Figure2.tif'), dpi=600, facecolor='white',
            pil_kwargs={'compression': 'tiff_lzw'})
print('Figure 2 written:', len(ORDER), 'species,', len(ROWS), 'comparisons')
