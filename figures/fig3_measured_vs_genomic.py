"""Figure 3 — the measured record set against the genomic record.

Panel A puts the reported D-alanine value and the layer-by-layer reporting
grid on one strain axis. Strains are grouped by predicted removal capability;
a third group holds those characterised without a reported value. Open
symbols mark strains with no deposited genome, which are excluded from the
test in Table S11d. One strain carries a value but is absent from the layer
compilation and is marked as such. Panel B sets genomes screened per species
against strains for which a value has ever been measured. Source data:
data/tables/TableS1_structure_map.csv and figures/fig2_data.json.
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

fig = plt.figure(figsize=(174 * MM, 182 * MM), dpi=600)

# ================= PANEL A : one strain axis, value and reporting grid =========
import csv as _csv

S1 = os.path.join(HERE, '..', 'data', 'tables', 'TableS1_structure_map.csv')
LAYERS = [('D-alanine (%)', 'D-alanine'), ('Chain length (n)', 'chain length'),
          ('Glycosyl (%)', 'glycosyl'), ('Anchor sugars (n)', 'anchor sugars'),
          ('Acyl chains (n)', 'acyl chains')]
NONLAB = {'Staphylococcus', 'Listeria', 'Streptococcus'}

# short name in the value series  ->  strain as written in Table S1
MAP = {
    'Lc. rhamnosus GG': 'Lacticaseibacillus rhamnosus GG (ATCC 53103)',
    'Lc. casei BL23': 'Lacticaseibacillus casei BL23',
    'Lp. plantarum L-137': 'L. plantarum L-137',
    'Lp. plantarum NCIMB 8826': 'L. plantarum NCIMB 8826',
    'Lc. rhamnosus DSM 20021': 'L. rhamnosus DSM 20021T',
    'Lp. pentosus DSM 20314': 'Lactiplantibacillus pentosus DSM 20314T',
    'Li. reuteri 100-23': 'Limosilactobacillus reuteri 100-23',
    'Lb. helveticus JCM 1120': 'Lactobacillus helveticus JCM 1120T',
    'Lb. gasseri JCM 1131': 'Lactobacillus gasseri JCM 1131T',
    'Lb. gasseri JCM 5814': 'Lactobacillus gasseri JCM 5814',
    'Lb. gasseri VLG3': 'Lactobacillus gasseri VLG3',
    'Lb. paragasseri JCM 1130': 'Lactobacillus paragasseri JCM 1130',
    'Lb. paragasseri VLG1': 'Lactobacillus paragasseri VLG1',
    'Lb. paragasseri VLG2': 'Lactobacillus paragasseri VLG2',
    'Lb. paragasseri VLG4': 'Lactobacillus paragasseri VLG4',
    'Lb. delbrueckii Ads-5': 'L. delbrueckii subsp. lactis Ads-5',
    'Lb. delbrueckii LL78': 'L. delbrueckii subsp. lactis LL78',
    'Lb. delbrueckii ATCC 15808': 'L. delbrueckii subsp. lactis ATCC 15808',
    'Lt. curvatus CP2998': 'Latilactobacillus curvatus CP2998',
    'Ap. kunkeei JCM 16173': 'Apilactobacillus kunkeei JCM 16173T',
}


def _state(v):
    v = (v or '').strip()
    if not v or v.upper() in ('NR', 'ND', 'NA', '-'):
        return 0
    return 2 if '-' in v else 1


def _short(name):
    parts = name.split()
    ab = {'Lactobacillus': 'Lb.', 'Lactiplantibacillus': 'Lp.',
          'Lacticaseibacillus': 'Lc.', 'Limosilactobacillus': 'Li.',
          'Latilactobacillus': 'Lt.', 'Levilactobacillus': 'Lv.',
          'Apilactobacillus': 'Ap.', 'L.': 'L.'}
    return (ab.get(parts[0], parts[0][:3] + '.') + ' ' + ' '.join(parts[1:])).replace('subsp. ', '')


_s1 = {r['Strain']: r for r in _csv.DictReader(open(S1, encoding='utf-8-sig'))
       if (r.get('Strain') or '').strip() and not r['Strain'].startswith('Note')
       and r['Strain'].split()[0] not in NONLAB}

has = [s for s in S if s[1] == 'HAS']
no = [s for s in S if s[1] == 'NO']
measured = {s[0] for s in S}
rest = [n for n in _s1 if n not in set(MAP.values())]
rest.sort(key=lambda n: -sum(1 for c, _ in LAYERS if _state(_s1[n].get(c))))

ROWS = ([('H', s) for s in has] + [None]
        + [('N', s) for s in no] + [None]
        + [('R', n) for n in rest])

n_rows = len([r for r in ROWS if r])
GAP = 0.9
H_ROW = 1.0

ax = fig.add_axes([0.255, 0.505, 0.300, 0.455])     # value plot
dx = fig.add_axes([0.590, 0.505, 0.190, 0.455])     # reporting grid

y, ticks, labels, bands = 0, [], [], {'H': [], 'N': [], 'R': []}
GRID = {0: '#F2F2F2', 1: NAVY, 2: '#9FB6D4'}
for item in ROWS:
    if item is None:
        y += GAP
        continue
    grp, r = item
    if grp in ('H', 'N'):
        name, cls, lo, hi, ext, dep = r
        c = NAVY if ext == 'BuOH' else RED
        if hi > lo:
            ax.plot([lo, hi], [y, y], color=c, lw=1.6, solid_capstyle='round', zorder=3)
        ax.plot([(lo + hi) / 2], [y], 'o' if dep else 's', ms=4.4, color=c,
                mfc=c if dep else 'white', mew=1.0, zorder=4)
        key = MAP.get(name)
        labels.append(name)
    else:
        key = r
        labels.append(_short(r))
    # reporting grid
    if key in _s1:
        for k, (col, _) in enumerate(LAYERS):
            dx.add_patch(Rectangle((k + 0.08, y - 0.38), 0.84, 0.76,
                                   fc=GRID[_state(_s1[key].get(col))], ec='white', lw=.5))
        nrep = sum(1 for col, _ in LAYERS if _state(_s1[key].get(col)))
        dx.text(5.32, y, str(nrep), ha='center', va='center', fontsize=5.4,
                color=NAVY if nrep == 5 else TXT, fontweight='bold' if nrep == 5 else 'normal')
    else:
        dx.text(2.5, y, 'not in the layer compilation', ha='center', va='center',
                fontsize=5.2, color=GREY, style='italic')
    ticks.append(y); bands[grp].append(y)
    y += H_ROW

TOP_Y = -0.8
BOT_Y = y - 0.4
for a_ in (ax, dx):
    a_.set_ylim(BOT_Y, TOP_Y)

ax.set_yticks(ticks)
ax.set_yticklabels(labels, fontsize=5.4, color=TXT, style='italic')
ax.set_xlim(-4, 88)
ax.set_xticks([0, 20, 40, 60, 80])
ax.set_xlabel('D-alanine substitution reported (%)', fontsize=6.8, color=TXT, labelpad=3)
ax.tick_params(labelsize=6.4, colors=TXT, length=3.0, width=1.0, color='#4D4D4D', pad=2)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_color('#4D4D4D'); ax.spines[sp].set_linewidth(1.0)
ax.grid(axis='x', color='#EDEDED', lw=.6, zorder=0); ax.set_axisbelow(True)

dx.set_xlim(0, 5.7)
dx.set_yticks([])
dx.set_xticks([k + 0.5 for k in range(5)])
dx.set_xticklabels([lab for _, lab in LAYERS], fontsize=5.4, color=TXT, rotation=38,
                   ha='right', rotation_mode='anchor')
dx.text(5.32, BOT_Y + 0.55, 'n', ha='center', va='center', fontsize=5.4,
        color=TXT, style='italic')
for sp in dx.spines.values():
    sp.set_visible(False)
dx.tick_params(length=0, pad=2)
dx.add_patch(Rectangle((1.02, TOP_Y + 0.35), 0.96, BOT_Y - TOP_Y - 0.7, fc='none',
                       ec=RED, lw=1.0, zorder=5))

# group bands down the left edge
for grp, col, lab in (('H', NAVY, 'removal predicted'),
                      ('N', RED, 'no removal predicted'),
                      ('R', GREY, 'no value reported')):
    ys = bands[grp]
    if not ys:
        continue
    y0, y1 = min(ys) - 0.5, max(ys) + 0.5
    ax.add_patch(Rectangle((-72, y0), 186, y1 - y0,
                           fc='#EEF2F8' if grp == 'H' else 'none', ec='none',
                           zorder=0, clip_on=False))
    ax.text(-70, (y0 + y1) / 2, lab, ha='center', va='center', rotation=90,
            fontsize=6.0, color=col, fontweight='bold', clip_on=False)

leg = [Line2D([], [], color=NAVY, lw=1.6, marker='o', ms=4.4, label='butanol extraction'),
       Line2D([], [], color=RED, lw=1.6, marker='o', ms=4.4, label='phenol extraction'),
       Line2D([], [], color=GREY, lw=0, marker='s', ms=4.4, mfc='white', mew=1.0,
              label='no deposited genome')]
lg = ax.legend(handles=leg, loc='lower left', bbox_to_anchor=(-0.02, 1.055), ncol=3,
               frameon=False, fontsize=5.6, handlelength=1.3, borderpad=0, columnspacing=0.9)
for t in lg.get_texts():
    t.set_color(TXT)

lgC = [Rectangle((0, 0), 1, 1, fc=GRID[1], ec='none'),
       Rectangle((0, 0), 1, 1, fc=GRID[2], ec='none'),
       Rectangle((0, 0), 1, 1, fc=GRID[0], ec='none')]
lg2 = dx.legend(lgC, ['point value', 'range only', 'not reported'],
                loc='lower right', bbox_to_anchor=(1.02, 1.005), ncol=3, frameon=False,
                fontsize=5.6, handlelength=0.9, handleheight=0.85, handletextpad=0.35,
                columnspacing=0.8, borderpad=0)
for t in lg2.get_texts():
    t.set_color(TXT)

fig.text(0.795, 0.505 + 0.455 * 0.66,
         f'{len(_s1)} strains have had their\nlipoteichoic acid characterised.\n'
         f'{sum(1 for n in _s1 if all(_state(_s1[n].get(c)) for c, _ in LAYERS))} report all five layers.\n'
         f'{sum(1 for n in _s1 if not _state(_s1[n].get("Chain length (n)")))} omit chain length \u2014 the layer\n'
         'Table 1 marks as changing potency\nand as unreadable from sequence.',
         ha='left', va='top', fontsize=5.8, color=TXT, linespacing=1.5)
fig.text(0.795, 0.505 + 0.455 * 0.30,
         'The values are not in one currency.\n'
         'The primary sources report it as mol\n'
         'per mol phosphate, as mol per mol\n'
         'glycerol, and as per cent of repeat\n'
         'units. The column here is the\ncompilation after conversion, not\nwhat was published.',
         ha='left', va='top', fontsize=5.4, color=RED, linespacing=1.5)

# ================= PANEL B =================
MEAS = {'Lb. gasseri': 3, 'Lb. paragasseri': 4, 'Lb. delbrueckii': 3, 'Lp. plantarum': 3,
        'Lc. rhamnosus': 2, 'Lb. helveticus': 1, 'Li. reuteri': 1, 'Lt. curvatus': 1,
        'Lc. casei': 1, 'Lp. pentosus': 1, 'Ap. kunkeei': 1}
order = sorted(D, key=lambda s: (-MEAS.get(s, 0), -D[s]['n']))
yy = np.arange(len(order))

bx = fig.add_axes([0.235, 0.075, 0.470, 0.330])
cx = fig.add_axes([0.735, 0.075, 0.105, 0.330])

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
fig.text(0.855, 0.075 + 0.330 * (1 - (n_meas + len(order)) / (2 * len(order))),
         f'{n_zero} of 30 species\n({g_zero:,} genomes)\nhave no measured\nvalue at all',
         ha='left', va='center', fontsize=6.2, color=RED, fontweight='bold', linespacing=1.5)

fig.text(0.020, 0.975, 'A', ha='left', va='top', fontsize=9, fontweight='bold', color='black')
fig.text(0.020, 0.440, 'B', ha='left', va='top', fontsize=9, fontweight='bold', color='black')


fig.savefig(os.path.join(HERE, 'Figure3.png'), dpi=600, facecolor='white')
fig.savefig(os.path.join(HERE, 'Figure3.tif'), dpi=600, facecolor='white',
            pil_kwargs={'compression': 'tiff_lzw'})
print('n_zero_species', n_zero, 'genomes', g_zero)
print('panel A strains', len(S))
