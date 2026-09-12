"""Figure 3 — the measured record set against the genomic record.

Panel A draws every reported D-alanine substitution value for a lactic acid
bacterium within the thirty species carried through the localization screen,
grouped by predicted removal capability. Open symbols mark strains with no
deposited genome; those strains are excluded from the test in Table S11d.
Panel B sets genomes screened per species against strains for which a value
has ever been measured. Source data: data/tables/TableS1_structure_map.csv
and figures/fig2_data.json.
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np, json

F = 'Liberation Sans'
plt.rcParams.update({'font.family': F, 'font.sans-serif': [F],
                     'mathtext.fontset': 'custom', 'mathtext.rm': F,
                     'mathtext.it': F + ':italic', 'mathtext.default': 'regular',
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
NAVY = '#1F3864'; RED = '#C00000'; BLUE = '#5B8FC9'; TXT = '#333333'; GREY = '#9A9A9A'
MM = 1 / 25.4

import os
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'fig2_data.json')))
ABBR = {'Lac.': 'Lactococcus', 'Li.': 'Limosilactobacillus', 'Lb.': 'Lactobacillus',
        'Lg.': 'Ligilactobacillus', 'Lt.': 'Latilactobacillus', 'P.': 'Pediococcus',
        'Lc.': 'Lacticaseibacillus', 'Lv.': 'Levilactobacillus', 'Lp.': 'Lactiplantibacillus',
        'Leu.': 'Leuconostoc', 'Ap.': 'Apilactobacillus'}
CLS = {}
for sp, v in D.items():
    CLS[sp] = 'HAS' if (v['counts'][0] + v['counts'][1]) / v['n'] >= 0.5 else 'NO'

# ---- panel A : every measured D-alanine value in the lactic acid bacteria ----
# (strain label, class, low, high, extraction, genome deposited)
S = [
 ('Lc. rhamnosus GG',            'HAS', 72, 74, 'BuOH',   True),
 ('Lc. casei BL23',              'HAS', 64, 64, 'BuOH',   True),
 ('Lp. plantarum L-137',         'HAS', 50, 50, 'BuOH',   False),
 ('Lp. plantarum ATCC 14917',    'HAS', 42, 42, 'BuOH',   True),
 ('Lp. plantarum NCIMB 8826',    'HAS', 42, 42, 'BuOH',   True),
 ('Lc. rhamnosus DSM 20021',     'HAS', 0,  0,  'phenol', True),
 ('Lp. pentosus DSM 20314',      'HAS', 0,  0,  'phenol', True),

 ('Li. reuteri 100-23',          'NO', 74, 79, 'BuOH',   True),
 ('Lb. helveticus JCM 1120',     'NO', 57, 64, 'phenol', True),
 ('Lb. gasseri JCM 1131',        'NO', 55, 55, 'BuOH',   True),
 ('Lb. gasseri JCM 5814',        'NO', 40, 55, 'BuOH',   False),
 ('Lb. gasseri VLG3',            'NO', 40, 55, 'BuOH',   True),
 ('Lb. paragasseri JCM 1130',    'NO', 40, 55, 'BuOH',   True),
 ('Lb. paragasseri VLG2',        'NO', 40, 55, 'BuOH',   True),
 ('Lb. delbrueckii Ads-5',       'NO', 42, 49, 'phenol', False),
 ('Lb. delbrueckii LL78',        'NO', 25, 28, 'phenol', False),
 ('Lb. paragasseri VLG4',        'NO', 24, 24, 'BuOH',   True),
 ('Lb. delbrueckii ATCC 15808',  'NO', 21, 27, 'phenol', False),
 ('Lt. curvatus CP2998',         'NO', 18, 18, 'BuOH',   False),
 ('Lb. paragasseri VLG1',        'NO', 17, 17, 'BuOH',   True),
 ('Ap. kunkeei JCM 16173',       'NO', 0,  0,  'BuOH',   True),
]

fig = plt.figure(figsize=(174 * MM, 228 * MM), dpi=600)

# ================= PANEL A =================
ax = fig.add_axes([0.235, 0.740, 0.720, 0.222])
has = [s for s in S if s[1] == 'HAS']
no  = [s for s in S if s[1] == 'NO']
rows = has + [None] + no
y = 0
ticks, labels = [], []
for r in rows:
    if r is None:
        y += 0.9
        continue
    name, cls, lo, hi, ext, dep = r
    c = NAVY if ext == 'BuOH' else RED
    if hi > lo:
        ax.plot([lo, hi], [y, y], color=c, lw=1.6, solid_capstyle='round', zorder=3)
    ax.plot([(lo + hi) / 2], [y], 'o' if dep else 's', ms=4.6, color=c,
            mfc=c if dep else 'white', mew=1.0, zorder=4)
    ticks.append(y); labels.append(name)
    y += 1

ax.set_ylim(y - 0.4, -0.8)
ax.set_xlim(-4, 88)
ax.set_yticks(ticks)
ax.set_yticklabels(labels, fontsize=5.6, color=TXT, style='italic')
ax.set_xticks([0, 20, 40, 60, 80])
ax.set_xlabel('D-alanine substitution reported (%)', fontsize=7.0, color=TXT, labelpad=3)
ax.tick_params(labelsize=6.6, colors=TXT, length=3.0, width=1.0, color='#4D4D4D', pad=2)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_color('#4D4D4D'); ax.spines[sp].set_linewidth(1.0)
ax.grid(axis='x', color='#EDEDED', lw=.6, zorder=0); ax.set_axisbelow(True)

nh = len(has)
ax.add_patch(Rectangle((-34, -0.55), 150, nh, fc='#EEF2F8', ec='none',
                       zorder=0, clip_on=False))
ax.text(-33, (nh - 1) / 2, 'removal predicted', ha='center', va='center', rotation=90,
        fontsize=6.2, color=NAVY, fontweight='bold', clip_on=False)
ax.text(-33, nh + 0.9 + (len(no) - 1) / 2, 'no removal predicted', ha='center', va='center',
        rotation=90, fontsize=6.2, color=RED, fontweight='bold', clip_on=False)

leg = [Line2D([], [], color=NAVY, lw=1.6, marker='o', ms=4.6, label='butanol extraction'),
       Line2D([], [], color=RED, lw=1.6, marker='o', ms=4.6, label='phenol extraction'),
       Line2D([], [], color=GREY, lw=0, marker='s', ms=4.6, mfc='white', mew=1.0,
              label='no deposited genome')]
lg = ax.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.5, 1.005), ncol=3,
               frameon=False, fontsize=6.2, handlelength=1.5, borderpad=0, columnspacing=1.6)
for t in lg.get_texts():
    t.set_color(TXT)

# ================= PANEL B =================
MEAS = {'Lb. gasseri': 3, 'Lb. paragasseri': 4, 'Lb. delbrueckii': 3, 'Lp. plantarum': 3,
        'Lc. rhamnosus': 2, 'Lb. helveticus': 1, 'Li. reuteri': 1, 'Lt. curvatus': 1,
        'Lc. casei': 1, 'Lp. pentosus': 1, 'Ap. kunkeei': 1}
order = sorted(D, key=lambda s: (-MEAS.get(s, 0), -D[s]['n']))
yy = np.arange(len(order))

bx = fig.add_axes([0.235, 0.408, 0.470, 0.268])
cx = fig.add_axes([0.735, 0.408, 0.105, 0.268])

bx.barh(yy, [D[s]['n'] for s in order], height=0.60, color='#D6DEEA', zorder=2)
bx.set_yticks(yy)
bx.set_yticklabels(order, fontsize=5.0, color=TXT, style='italic')
bx.set_ylim(len(order) - 0.45, -0.75)
bx.set_xlim(0, 520)
bx.set_xticks([0, 250, 500])
bx.set_xlabel('Genomes screened', fontsize=6.8, color=TXT, labelpad=3)

cx.barh(yy, [MEAS.get(s, 0) for s in order], height=0.60,
        color=[NAVY if CLS[s] == 'HAS' else RED for s in order], zorder=2)
cx.set_yticks(yy); cx.set_yticklabels([])
cx.set_ylim(len(order) - 0.45, -0.75)
cx.set_xlim(0, 4.6); cx.set_xticks([0, 2, 4])
cx.set_xlabel('Strains with a\nmeasured value', fontsize=6.8, color=TXT, labelpad=3,
              linespacing=1.3)

for a in (bx, cx):
    a.tick_params(labelsize=6.2, colors=TXT, length=3.0, width=1.0, color='#4D4D4D', pad=2)
    for sp in ('top', 'right'):
        a.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        a.spines[sp].set_color('#4D4D4D'); a.spines[sp].set_linewidth(1.0)
    a.grid(axis='x', color='#EDEDED', lw=.6, zorder=0); a.set_axisbelow(True)

n_meas = len(MEAS)
cx.plot([-0.35, -0.35], [-0.45, n_meas - 0.55], color=TXT, lw=1.0, clip_on=False)
cx.plot([-0.35, -0.35], [n_meas - 0.45, len(order) - 0.55], color=GREY, lw=1.0, clip_on=False)
n_zero = len(order) - n_meas
g_zero = sum(D[s]['n'] for s in order if s not in MEAS)
cx.add_patch(Rectangle((0, n_meas - 0.45), 4.6, len(order) - n_meas,
                       fc='#F6EBEB', ec='none', zorder=0))
fig.text(0.855, 0.408 + 0.268 * (1 - (n_meas + len(order)) / (2 * len(order))),
         f'{n_zero} of 30 species\n({g_zero:,} genomes)\nhave no measured\nvalue at all',
         ha='left', va='center', fontsize=6.2, color=RED, fontweight='bold', linespacing=1.5)

fig.text(0.020, 0.975, 'A', ha='left', va='top', fontsize=9, fontweight='bold', color='black')
fig.text(0.020, 0.500, 'B', ha='left', va='top', fontsize=9, fontweight='bold', color='black')


# ================= PANEL C : what has actually been reported =================
import csv as _csv, re as _re
S1 = os.path.join(HERE, '..', 'data', 'tables', 'TableS1_structure_map.csv')
LAYERS = [('D-alanine (%)', 'D-alanine'), ('Chain length (n)', 'chain length'),
          ('Glycosyl (%)', 'glycosyl'), ('Anchor sugars (n)', 'anchor sugars'),
          ('Acyl chains (n)', 'acyl chains')]
NONLAB = {'Staphylococcus', 'Listeria', 'Streptococcus'}


def _state(v):
    v = (v or '').strip()
    if not v or v.upper() in ('NR', 'ND', 'NA', '-'):
        return 0
    return 2 if '-' in v else 1


def _short(name):
    parts = name.split()
    g, rest = parts[0], ' '.join(parts[1:])
    ab = {'Lactobacillus': 'Lb.', 'Lactiplantibacillus': 'Lp.', 'Lacticaseibacillus': 'Lc.',
          'Limosilactobacillus': 'Li.', 'Latilactobacillus': 'Lt.', 'Levilactobacillus': 'Lv.',
          'Apilactobacillus': 'Ap.', 'L.': 'L.'}
    return f"{ab.get(g, g[:3] + '.')} {rest}"


_rows = [r for r in _csv.DictReader(open(S1, encoding='utf-8-sig'))
         if (r.get('Strain') or '').strip() and not r['Strain'].startswith('Note')
         and r['Strain'].split()[0] not in NONLAB]
_rows.sort(key=lambda r: (-sum(1 for c, _ in LAYERS if _state(r.get(c))), r['Strain']))

dx = fig.add_axes([0.300, 0.052, 0.300, 0.278])
GRID = {0: '#F2F2F2', 1: NAVY, 2: '#9FB6D4'}
for i, r in enumerate(_rows):
    for j, (col, _) in enumerate(LAYERS):
        dx.add_patch(Rectangle((j + 0.08, i + 0.12), 0.84, 0.76,
                               fc=GRID[_state(r.get(col))], ec='white', lw=.5))
    n_rep = sum(1 for c, _ in LAYERS if _state(r.get(c)))
    dx.text(5.28, i + 0.5, str(n_rep), ha='center', va='center', fontsize=5.4,
            color=NAVY if n_rep == 5 else TXT, fontweight='bold' if n_rep == 5 else 'normal')
dx.set_xlim(0, 5.6)
dx.set_ylim(len(_rows), 0)
dx.set_yticks([i + 0.5 for i in range(len(_rows))])
dx.set_yticklabels([_short(r['Strain'])[:30] for r in _rows], fontsize=4.8,
                   color=TXT, style='italic')
dx.set_xticks([j + 0.5 for j in range(5)])
dx.set_xticklabels([lab for _, lab in LAYERS], fontsize=5.4, color=TXT, rotation=38,
                   ha='right', rotation_mode='anchor')
dx.text(5.28, -0.55, 'n', ha='center', va='center', fontsize=5.4, color=TXT, style='italic')
for sp in dx.spines.values():
    sp.set_visible(False)
dx.tick_params(length=0, pad=2)

n_full = sum(1 for r in _rows if all(_state(r.get(c)) for c, _ in LAYERS))
n_nochain = sum(1 for r in _rows if not _state(r.get('Chain length (n)')))
dx.add_patch(Rectangle((1.02, -0.05), 0.96, len(_rows) + 0.1, fc='none', ec=RED,
                       lw=1.0, zorder=5))
fig.text(0.617, 0.300,
         f'{len(_rows)} strains have had their\nlipoteichoic acid characterised.\n'
         f'{n_full} report all five layers.\n{n_nochain} omit chain length \u2014 the layer\n'
         'Table 1 marks as changing potency\nand as unreadable from sequence.',
         ha='left', va='top', fontsize=6.0, color=TXT, linespacing=1.5)

lgC = [Rectangle((0, 0), 1, 1, fc=GRID[1], ec='none'),
       Rectangle((0, 0), 1, 1, fc=GRID[2], ec='none'),
       Rectangle((0, 0), 1, 1, fc=GRID[0], ec='none')]
lg3 = dx.legend(lgC, ['point value', 'range only', 'not reported'],
                loc='lower center', bbox_to_anchor=(0.5, 1.035), ncol=3, frameon=False,
                fontsize=5.8, handlelength=1.0, handleheight=0.85, handletextpad=0.42,
                columnspacing=1.2, borderpad=0)
for t in lg3.get_texts():
    t.set_color(TXT)

fig.text(0.617, 0.140,
         'And the values are not in one currency. The primary\n'
         'sources report substitution three ways: mol D-alanine\n'
         'per mol phosphate (MacArthur & Archibald 1984),\n'
         'mol alanyl-glycerol per mol total glycerol (Fischer &\n'
         'R\u00f6sel 1980), and per cent of repeat units substituted\n'
         '(Shiraishi et al. 2025). The D-alanine column here is\n'
         'the compilation after conversion, not what was published.',
         ha='left', va='top', fontsize=5.8, color=RED, linespacing=1.5)

fig.text(0.020, 0.352, 'C', ha='left', va='top', fontsize=9, fontweight='bold', color='black')

fig.savefig(os.path.join(HERE, 'Figure3.png'), dpi=600, facecolor='white')
fig.savefig(os.path.join(HERE, 'Figure3.tif'), dpi=600, facecolor='white',
            pil_kwargs={'compression': 'tiff_lzw'})
print('n_zero_species', n_zero, 'genomes', g_zero)
print('panel A strains', len(S))
