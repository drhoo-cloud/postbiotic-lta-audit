import os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, csv
F='Liberation Sans'
plt.rcParams.update({'font.family':F,'font.sans-serif':[F],'mathtext.fontset':'custom',
 'mathtext.rm':F,'mathtext.it':F+':italic','mathtext.default':'regular',
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
NAVY='#1F3864'; RED='#C00000'; TXT='#333333'; AX='#9A9A9A'
ACIDO={'Lactobacillus crispatus','Lactobacillus helveticus','Lactobacillus gasseri','Lactobacillus jensenii',
 'Lactobacillus acidophilus','Lactobacillus johnsonii','Lactobacillus paragasseri','Lactobacillus amylovorus'}
sc=[]; ac=[]
SRC=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','tables','TableS9a.csv')
if not os.path.exists(SRC):
    raise SystemExit('figS1 needs data/tables/TableS9a.csv. Run scripts/05_aggregate.py first, '
                     'or download TableS9a from the Zenodo data record.')
for r in csv.DictReader(open(SRC,encoding='utf-8-sig')):
    if r['species'].split()[0] in ('Enterococcus','Bifidobacterium'): continue
    try: v=float(r['DltE_score'])
    except: continue
    if v<=0: continue
    sc.append(v)
    if r['species'] in ACIDO: ac.append(v)
sc=np.array(sc); ac=np.array(ac)
MM=1/25.4
fig=plt.figure(figsize=(114*MM,68*MM),dpi=600)
ax=fig.add_axes([0.135,0.185,0.845,0.735])
b=np.arange(105,325,3)
ax.hist(sc,bins=b,color='#C9D3E4',ec='white',lw=.3,zorder=3,label='all lactic acid bacteria')
ax.hist(ac,bins=b,color=RED,ec='white',lw=.3,zorder=4,label='$\\it{L.\\ acidophilus}$ group')
for t,st in ((110,'-'),(130,(0,(4,2)))):
    ax.axvline(t,color=NAVY,lw=1.0 if t==110 else .8,ls=st,zorder=5)
    ax.text(t,ax.get_ylim()[1],'',ha='center')
ax.set_xlim(105,322)
ax.set_xlabel('DltE alignment score (HMMER bit score)',fontsize=7.0,color=TXT,labelpad=2)
ax.set_ylabel('Genomes',fontsize=7.0,color=TXT,labelpad=3)
ax.tick_params(labelsize=6.8,colors=TXT,length=2.5,pad=2)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(AX); ax.spines[s].set_linewidth(.7)
ax.grid(axis='y',color='#EDEDED',lw=.6,zorder=0); ax.set_axisbelow(True)
ymax=ax.get_ylim()[1]
for t,lab,yy in ((110,'110  threshold used',0.985),(130,'130  sensitivity check',0.885)):
    ax.text(t+2.5,ymax*yy,lab,ha='left',va='top',fontsize=6.4,color=NAVY,fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.14',fc='white',ec='none'))
lo,hi=ac.min(),ac.max()
ax.annotate('divergent 308-314 aa form of the\nL. acidophilus group, scores %.0f-%.0f'%(lo,hi),xy=(128,ymax*0.36),
            xytext=(196,ymax*0.66),fontsize=6.6,color=RED,ha='left',va='center',linespacing=1.35,
            arrowprops=dict(arrowstyle='->',color=RED,lw=.7))
leg=ax.legend(loc='upper right',fontsize=6.4,frameon=False,handlelength=1.0,handleheight=.8,borderpad=0)
for t in leg.get_texts(): t.set_color(TXT)
fig.text(0.985,0.035,'686 genomes score above 320 and are not shown; Table S9a records no score below 110',ha='right',va='bottom',fontsize=5.8,color='#8A8A8A')
_H=os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(_H,'FigureS1.png'),dpi=600,facecolor='white')
fig.savefig(os.path.join(_H,'FigureS1.tif'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
print('n=%d  acido n=%d  range %.1f-%.1f'%(len(sc),len(ac),lo,hi))
