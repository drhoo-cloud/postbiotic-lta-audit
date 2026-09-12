import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
F='Liberation Sans'
plt.rcParams.update({'font.family':F,'font.sans-serif':[F],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
NAVY='#1F3864'; RED='#C00000'; TXT='#333333'; AX='#9A9A9A'
RES=['Ser128','Tyr213','Tyr338','Gly348','His347','Arg345','Leu349','Leu127']
CONS=[100.0,100.0,99.9,99.6,39.2,19.4,11.4,8.1]
DIST=['2.4 \u00c5','2.8 \u00c5','2.9 \u00c5','2.5 \u00c5','3.2 \u00c5','3.8 \u00c5','2.0 \u00c5','3.0 \u00c5']
MM=1/25.4
fig=plt.figure(figsize=(72*MM,56*MM),dpi=600)
ax=fig.add_axes([0.190,0.300,0.780,0.590])
x=np.arange(8); cols=[NAVY]*4+[RED]*4
ax.bar(x,CONS,width=0.66,color=cols,zorder=3)
for xi,v,c in zip(x,CONS,cols):
    if v>60: ax.text(xi,v-9,f'{v:.0f}%',ha='center',va='top',fontsize=5.8,color='white',fontweight='bold')
    else:    ax.text(xi,v+3.5,f'{v:.1f}%',ha='center',va='bottom',fontsize=5.8,color=c,fontweight='bold')
ax.set_ylim(0,150); ax.set_xlim(-0.72,7.72); ax.set_yticks([0,25,50,75,100])
ax.set_ylabel('Conservation (%)',fontsize=6.4,color=TXT,labelpad=3)
ax.set_xticks(x); ax.set_xticklabels([f'{r}\n{d}' for r,d in zip(RES,DIST)],fontsize=5.2,color=TXT,linespacing=1.3)
ax.tick_params(labelsize=5.8,colors=TXT,length=2.5,pad=2)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'):
    ax.spines[s].set_color('#4D4D4D'); ax.spines[s].set_linewidth(1.1)
ax.tick_params(axis='both',length=3.5,width=1.0,color='#4D4D4D')
ax.grid(axis='y',color='#EDEDED',lw=.6,zorder=0); ax.set_axisbelow(True)
ax.axvline(3.5,color=AX,lw=.7,ls=(0,(4,3)),zorder=2)
ax.text(1.5,140,'catalytic / scaffold\nmean 99.9%',ha='center',va='center',fontsize=6.0,color=NAVY,fontweight='bold',linespacing=1.3)
ax.text(5.5,140,'substrate contact\nmean 19.5%',ha='center',va='center',fontsize=6.0,color=RED,fontweight='bold',linespacing=1.3)
ax.text(3.5,-38,'727 unique DltE proteins',ha='center',va='top',fontsize=5.8,color='#8A8A8A')
import os as _os
_H=_os.path.dirname(_os.path.abspath(__file__))
fig.savefig(_os.path.join(_H,'FigureS2.png'),dpi=600,facecolor='white')
fig.savefig(_os.path.join(_H,'FigureS2.tif'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})

print('ok')
