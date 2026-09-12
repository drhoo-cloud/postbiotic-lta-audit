import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
F='Liberation Sans'
plt.rcParams.update({'font.family':F,'font.sans-serif':[F],'mathtext.fontset':'custom','mathtext.rm':F,
 'mathtext.it':F+':italic','mathtext.default':'regular','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
PURPLE='#7B5EA7'; RED='#C00000'; TXT='#333333'; AX='#4D4D4D'
SP=['Lb. crispatus','Lb. helveticus','Lb. gasseri','Lb. jensenii','Lb. acidophilus','Lb. johnsonii',
 'Lb. paragasseri','Lb. amylovorus','Lac. lactis','Lac. garvieae','Lb. delbrueckii','P. acidilactici',
 'P. pentosaceus','Leu. mesenteroides','Lp. pentosus','Lc. rhamnosus','Lv. brevis','Li. fermentum',
 'Li. mucosae','Lg. salivarius','Lt. curvatus','Lp. plantarum','Lc. paracasei','Lt. sakei','Li. reuteri',
 'Leu. citreum','Lg. animalis','Ap. kunkeei','Lc. casei','Lp. paraplantarum']
LONG=[0]*13+[1]*17
SHORT=[1]*30
for i in (8,9): SHORT[i]=0
MM=1/25.4
fig=plt.figure(figsize=(100*MM,56*MM),dpi=600)
ax=fig.add_axes([0.205,0.330,0.775,0.470])
ax.set_xlim(-0.6,29.6); ax.set_ylim(2.25,-0.70); ax.axis('off')
SX,SY=0.68,0.44
for i in range(30):
    for row,vec in ((0,LONG),(1,SHORT)):
        on=vec[i]
        ax.add_patch(Rectangle((i-SX/2,row-SY/2),SX,SY,fc=PURPLE if on else 'white',
                               ec=PURPLE if on else '#C9C9C9',lw=.5,zorder=3))
    ax.text(i,1.94,SP[i],rotation=90,ha='center',va='top',fontsize=4.5,color=TXT,style='italic')
# axis lines
ax.plot([-0.6,29.6],[1.30,1.30],color=AX,lw=1.0,clip_on=False)
ax.plot([-0.62,-0.62],[-0.28,1.30],color=AX,lw=1.0,clip_on=False)
ax.text(-1.3,0,'long form',ha='right',va='center',fontsize=6.2,color=PURPLE,fontweight='bold')
ax.text(-1.3,0.30,'sets chain length',ha='right',va='center',fontsize=5.2,color=TXT)
ax.text(-1.3,1.0,'short form',ha='right',va='center',fontsize=6.2,color=TXT,fontweight='bold')
ax.text(-1.3,1.30,'in existing panels',ha='right',va='center',fontsize=5.2,color=TXT)
ax.plot([-0.55,12.55],[-0.50,-0.50],color=RED,lw=1.0)
ax.text(6.0,-0.64,'long form absent \u2014 13 of 30 species',ha='center',va='center',fontsize=6.0,color=RED,fontweight='bold')
fig.text(0.022,0.60,'MprF',ha='center',va='center',fontsize=7.0,color=TXT,fontweight='bold',rotation=90)
import os as _os
_H=_os.path.dirname(_os.path.abspath(__file__))
fig.savefig(_os.path.join(_H,'FigureS3.png'),dpi=600,facecolor='white')
fig.savefig(_os.path.join(_H,'FigureS3.tif'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
print('ok')
