import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
F='Liberation Sans'
plt.rcParams.update({'font.family':F,'font.sans-serif':[F],'mathtext.fontset':'custom','mathtext.rm':F,
 'mathtext.it':F+':italic','mathtext.default':'regular','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
NAVY='#1F3864'; RED='#C00000'; TXT='#333333'; BAR='#A8C4E0'
S=[('L. paragasseri','VLG1',17.0),('L. paragasseri','VLG4',24.0),('L. paragasseri','VLG2',47.5),
   ('L. gasseri','VLG3',47.5),('L. gasseri','JCM 1131',55.0)]
MM=1/25.4
fig=plt.figure(figsize=(140*MM,50*MM),dpi=600)
ax=fig.add_axes([0.150,0.255,0.800,0.660])
x=np.arange(5); v=[s[2] for s in S]
ax.vlines(x,0,v,color=BAR,lw=5,zorder=2)
ax.plot(x,v,'o',ms=6,color=NAVY,zorder=3)
for xi,vv in zip(x,v):
    ax.text(xi,vv+3.2,f'{vv}',ha='center',va='bottom',fontsize=7.6,color=NAVY,fontweight='bold')
ax.set_ylim(0,66); ax.set_xlim(-0.55,4.55)
ax.set_yticks([0,20,40,60])
ax.set_ylabel('D-alanine substitution (%)',fontsize=7.4,color=TXT,labelpad=4)
ax.set_xticks(x)
lab=[r'$\it{%s}$'%g.replace(' ','\\ ')+'\n'+st for g,st,_ in S]
ax.set_xticklabels(lab,fontsize=6.9,color=TXT,linespacing=1.6)
ax.tick_params(labelsize=7.0,colors=TXT,length=3.5,width=1.0,color='#4D4D4D',pad=3)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
for sp in ('left','bottom'): ax.spines[sp].set_color('#4D4D4D'); ax.spines[sp].set_linewidth(1.1)
ax.grid(axis='y',color='#EDEDED',lw=.6,zorder=0); ax.set_axisbelow(True)
ax.plot([0,4.35],[17,17],color=RED,lw=.8,ls=(0,(2,2)),zorder=2)
ax.plot([4,4.35],[55,55],color=RED,lw=.8,ls=(0,(2,2)),zorder=2)
ax.annotate('',xy=(4.42,55),xytext=(4.42,17),arrowprops=dict(arrowstyle='<->',color=RED,lw=1.3))
ax.text(4.30,36,'3.2-fold',ha='right',va='center',fontsize=7.6,color=RED,fontweight='bold')
ax.text(-0.45,62.5,'identical DltE and dlt proteins across all five',fontsize=6.8,color=TXT,style='italic')
fig.savefig('F1B_strains.png',dpi=600,facecolor='white')
print('ok')
