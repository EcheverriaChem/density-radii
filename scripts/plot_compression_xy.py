"""Wigner-Seitz radius against r_outer (Section 5.1.2, Figure 7).

Usage (from the package root):  python scripts/plot_compression_xy.py
"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':28,'axes.linewidth':2,
  'xtick.major.width':2,'ytick.major.width':2,'xtick.direction':'in','ytick.direction':'in',
  'xtick.top':True,'ytick.right':True,'lines.linewidth':2})
from pathlib import Path
D=Path(__file__).resolve().parent.parent/'results'/'revision'/'compression'
c=pd.read_csv(D/'compression_points.csv')
sets=[('Ambient',c.P_GPa==0,'o','#0060ad',True),
      ('s$\\rightarrow$d onset',c.state.str.contains('s-d'),'s','#dd181f',False),
      ('Electride',c.state.str.contains('electride'),'^','#00a000',True)]
fig,ax=plt.subplots(figsize=(8,8)); xx=np.linspace(1.25,2.7,10)
for lab,m,mk,col,fill in sets:
    g=c[m]; x,y=g.r_outer.values,g.r_WS.values
    k=(x@y)/(x@x)                         # fit through the origin
    rms=np.sqrt(np.mean((y-k*x)**2))
    if len(x)>2:
        b=np.polyfit(x,y,1); r2=np.corrcoef(x,y)[0,1]**2
        print(f"{lab:22s} n={len(x)} k={k:.3f} rms={rms:.3f} A | free fit: slope={b[0]:.3f} intercept={b[1]:+.3f} r2={r2:.3f}")
    else: print(f"{lab:22s} n={len(x)} k={k:.3f} rms={rms:.3f} A")
    ax.plot(xx,k*xx,'--',color=col)
    ax.plot(x,y,mk,ms=14,mfc=col if fill else 'none',mec=col,mew=2,ls='none',label=lab)
    for e,xi,yi in zip(g.el,x,y): ax.annotate(e,(xi,yi),xytext=(8,-22),textcoords='offset points',fontsize=16)
ax.plot(xx,xx,':',color='grey')
ax.set_xlim(1.25,2.7); ax.set_ylim(0.7,3.1); ax.set_xticks([1.4,1.8,2.2,2.6])
ax.set_xlabel('$r_{\\rm outer}$ (\\u00c5)'.replace('\\\\u00c5','Å')); ax.set_ylabel('$r_{\\rm WS}$ (Å)')
ax.set_xlabel('$r_{\\rm outer}$ (Å)')
ax.legend(frameon=False,fontsize=22,loc='upper left')
fig.tight_layout(); fig.savefig(D/'Fig_compression_xy.png',dpi=150); fig.savefig(D/'Fig_compression_xy.eps')
