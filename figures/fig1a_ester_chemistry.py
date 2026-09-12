import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Ellipse, Rectangle, FancyBboxPatch
import numpy as np
F='Liberation Sans'
plt.rcParams.update({'font.family':F,'font.sans-serif':[F],'mathtext.fontset':'custom',
 'mathtext.rm':F,'mathtext.it':F+':italic','mathtext.bf':F+':bold','mathtext.default':'regular',
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
NAVY='#1F3864'; RED='#C00000'; BLUE='#5B8FC9'; GREEN='#2E7D5B'; TXT='#333333'; GREY='#9A9A9A'
POS='#C00000'; NEG='#1F3864'
MM=1/25.4
fig=plt.figure(figsize=(140*MM,74*MM),dpi=600)

# ================= PANEL A : the ester bond =================
A=fig.add_axes([0.0,0.02,1.0,0.96]); A.set_xlim(0,100); A.set_ylim(0,100); A.axis('off')
def bond(x1,y1,x2,y2,c=TXT,lw=1.0,ls='-'): A.plot([x1,x2],[y1,y2],color=c,lw=lw,ls=ls,zorder=2,solid_capstyle='round')
def atom(x,y,s,c=TXT,fs=7.4,w='normal'):
    A.text(x,y,s,ha='center',va='center',fontsize=fs,color=c,fontweight=w,zorder=4,
           bbox=dict(boxstyle='round,pad=0.12',fc='white',ec='none'))

def backbone(x0,y0,sub):
    """one glycerophosphate repeat; returns C2 coordinates"""
    xs=[x0,x0+5,x0+10,x0+15,x0+20,x0+25]
    ys=[y0,y0+4,y0,y0+4,y0,y0+4]
    for i in range(5): bond(xs[i],ys[i],xs[i+1],ys[i+1])
    atom(xs[0],ys[0],'O'); atom(xs[1],ys[1],'CH$_2$'); atom(xs[2],ys[2],'CH')
    atom(xs[3],ys[3],'CH$_2$'); atom(xs[4],ys[4],'O'); atom(xs[5],ys[5],'P',NAVY,8.0,'bold')
    # phosphate
    bond(xs[5],ys[5]+1.2,xs[5],ys[5]+7); bond(xs[5]-0.8,ys[5]+1.2,xs[5]-0.8,ys[5]+7)
    atom(xs[5],ys[5]+8.5,'O')
    bond(xs[5],ys[5]-1.2,xs[5],ys[5]-7); atom(xs[5],ys[5]-8.5,'O$^-$',NEG,7.6,'bold')
    bond(xs[5]+1.2,ys[5],xs[5]+6,ys[5]); atom(xs[5]+7.5,ys[5],'O')
    A.text(xs[5]+11,ys[5],'$n$',ha='center',va='center',fontsize=7.4,color=GREY,style='italic')
    bond(xs[0]-1.2,ys[0],xs[0]-5,ys[0],GREY,ls=(0,(2,2)))
    return xs[2],ys[2]

# --- left: free hydroxyl
A.text(8,92,'Unsubstituted backbone',ha='left',va='center',fontsize=7.6,color=NAVY,fontweight='bold')
cx,cy=backbone(11,66,False)
bond(cx,cy-1.4,cx,cy-8); atom(cx,cy-9.6,'OH',GREEN,7.6,'bold')
A.annotate('free hydroxyl\nat C2',xy=(cx+2.4,cy-9.6),xytext=(cx+9,cy-15.5),fontsize=6.8,color=GREEN,
           ha='left',va='center',linespacing=1.35,
           arrowprops=dict(arrowstyle='-',color=GREEN,lw=.6,shrinkA=1,shrinkB=1))

# --- D-alanine
A.text(8,38,'D-alanine',ha='left',va='center',fontsize=7.6,color=RED,fontweight='bold')
bx,by=15,26
bond(bx,by,bx+6,by+4); bond(bx+6,by+4,bx+12,by); bond(bx+6,by+4,bx+6,by+11)
atom(bx,by,'H$_3$N$^+$',POS,7.6,'bold'); atom(bx+6,by+4,'CH'); atom(bx+6,by+11,'CH$_3$')
atom(bx+12,by,'C',TXT,7.6)
bond(bx+12+1.0,by+0.6,bx+17,by+4); bond(bx+12+1.8,by-0.2,bx+17.8,by+3.2)
atom(bx+18.5,by+4.6,'O')
bond(bx+12+1.2,by-0.6,bx+17,by-4); atom(bx+18.5,by-5,'OH',GREEN,7.6,'bold')

# --- arrow
# --- who adds and who removes: the point of the paper, on the molecule
from matplotlib.patches import FancyArrowPatch as _FA
A.add_patch(_FA((40,50),(56,50),arrowstyle='-|>',mutation_scale=13,color=NAVY,lw=2.0))
A.text(48,60,'DltA\u2013D',ha='center',va='center',fontsize=8.4,color=NAVY,fontweight='bold')
A.text(48,54,'+ ATP',ha='center',va='center',fontsize=7.2,color=NAVY)
A.text(48,42,'four genes \u00b7 in every panel',ha='center',va='center',fontsize=6.9,color=NAVY)
A.add_patch(_FA((56,20),(40,20),arrowstyle='-|>',mutation_scale=13,color=RED,lw=2.0,linestyle=(0,(4,2))))
A.text(48,30,'DltE',ha='center',va='center',fontsize=8.4,color=RED,fontweight='bold')
A.text(48,24,'+ H$_2$O',ha='center',va='center',fontsize=7.2,color=RED)
A.text(48,12,'one gene \u00b7 in none',ha='center',va='center',fontsize=6.9,color=RED)


# --- right: D-alanyl ester
A.text(59,92,'D-alanyl ester',ha='left',va='center',fontsize=7.6,color=RED,fontweight='bold')
cx2,cy2=backbone(62,66,True)
bond(cx2,cy2-1.4,cx2,cy2-7); atom(cx2,cy2-8.4,'O',RED,7.6,'bold')
bond(cx2,cy2-9.8,cx2,cy2-15); atom(cx2,cy2-16.4,'C')
bond(cx2+1.0,cy2-15.8,cx2+5,cy2-13.4); bond(cx2+1.6,cy2-17.0,cx2+5.6,cy2-14.6)
atom(cx2+6.6,cy2-13.0,'O')
bond(cx2,cy2-17.8,cx2,cy2-23); atom(cx2,cy2-24.4,'CH')
bond(cx2-1.4,cy2-24.4,cx2-7,cy2-24.4); atom(cx2-9.0,cy2-24.4,'CH$_3$')
bond(cx2,cy2-25.8,cx2,cy2-31); atom(cx2,cy2-32.4,'H$_3$N$^+$',POS,7.6,'bold')
A.add_patch(Ellipse((cx2+1.6,cy2-12.2),13.5,10.5,angle=-18,fc='none',ec=RED,lw=.9,ls=(0,(3,2)),zorder=1))
A.text(cx2-9.5,cy2-12.0,'ester',ha='right',va='center',fontsize=7.0,color=RED,fontweight='bold')
A.text(97,cy2-30,'one $+$ added\nper substitution',ha='right',va='center',fontsize=6.8,color=POS,linespacing=1.35)


fig.savefig('F1A_chem.png',dpi=600,facecolor='white')

print('ok')
