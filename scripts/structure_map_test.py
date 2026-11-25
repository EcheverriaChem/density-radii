import numpy as np, pandas as pd, re, warnings
warnings.filterwarnings('ignore')
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import ShuffleSplit, LeaveOneOut
B2A=0.529177
O=pd.read_csv('octet_ghiringhelli2015.csv')
Zsym={3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',19:'K',20:'Ca',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',35:'Br',37:'Rb',38:'Sr',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',53:'I',55:'Cs',56:'Ba'}
O['A']=O.ZA.map(Zsym); O['B']=O.ZB.map(Zsym)
# our radii: valence s and p; empty valence p -> lowest virtual p
R=pd.read_csv('orbital_radii_by_l.csv'); R=R[R.charge==0].set_index('symbol')
def has_valence_p(conf):
    sub=conf.split(']')[-1].split()
    n=max([int(s[0]) for s in sub if s[1]=='s'] or [0])
    return any(s[1]=='p' and int(s[0])==n for s in sub)
ours=pd.DataFrame({'rs':R.r_s,'rp':[R.loc[s,'r_p'] if has_valence_p(R.loc[s,'configuration']) else R.loc[s,'rv_p'] for s in R.index]},index=R.index)
Zu=pd.read_csv('zunger1980_radii.csv').set_index('symbol'); zun=pd.DataFrame({'rs':Zu.rs_au*B2A,'rp':Zu.rp_au*B2A})
I=pd.read_csv('/home/claude/release/density-radii/results/ion_radii_selected.csv'); I=I[(I.method=='PBE0')&(I.scheme=='anorcc')]
IP=I[I.charge==1].set_index('symbol').IE_calc_eV; EA=I[I.charge==-1].set_index('symbol').EA_calc_eV.copy(); EA['N']=-0.07
y=O.dE_RS_ZB_eV.values; lab=(y>0).astype(int)
def feats(rad):
    a,b=rad.loc[O.A],rad.loc[O.B]
    rsig=np.abs((a.rp.values+a.rs.values)-(b.rp.values+b.rs.values)); rpi=np.abs(a.rp.values-a.rs.values)+np.abs(b.rp.values-b.rs.values)
    d1=(O.B.map(EA).values-O.B.map(IP).values)/a.rp.values**2; d2=np.abs(a.rs.values-b.rp.values)/np.exp(a.rs.values)
    return np.c_[rsig,rpi], np.c_[d1,d2], np.c_[a.rs,a.rp,b.rs,b.rp]
def lsq(X):
    A=np.c_[X,np.ones(len(X))]; c,*_=np.linalg.lstsq(A,y,rcond=None); e=A@c-y; return np.sqrt((e**2).mean()), np.abs(e).max()
def lsq_cv(X,n=150,seed=0):
    ss=ShuffleSplit(n_splits=n,test_size=0.1,random_state=seed); err=[]
    for tr,te in ss.split(X):
        A=np.c_[X,np.ones(len(X))]; c,*_=np.linalg.lstsq(A[tr],y[tr],rcond=None); err+=list(A[te]@c-y[te])
    return np.sqrt(np.mean(np.square(err)))
def krr_cv(X,n=150):
    Xs=(X-X.mean(0))/X.std(0); best=None
    for lam in [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]:
        for g in [0.01,0.03,0.1,0.3,1,3,10,30,100]:
            ss=ShuffleSplit(n_splits=n,test_size=0.1,random_state=1); err=[]
            for tr,te in ss.split(Xs):
                m=KernelRidge(alpha=lam,kernel='rbf',gamma=g).fit(Xs[tr],y[tr]); err+=list(m.predict(Xs[te])-y[te])
            r=np.sqrt(np.mean(np.square(err)))
            if best is None or r<best[0]: best=(r,lam,g)
    return best
def sep(X):
    Xs=(X-X.mean(0))/X.std(0)
    acc=LogisticRegression(C=1e4,max_iter=10000).fit(Xs,lab).score(Xs,lab)
    loo=np.mean([SVC(kernel='linear',C=100).fit(Xs[tr],lab[tr]).predict(Xs[te])[0]==lab[te][0] for tr,te in LeaveOneOut().split(Xs)])
    return acc,loo
print(f"82 compounds; {sum(lab==0)} RS (dE<0), {sum(lab==1)} ZB\n")
for name,rad in (('Zunger 1980 (pseudopotential)',zun),('this work (PBE0/DKH2 all-electron)',ours)):
    Z,G,R4=feats(rad)
    print(name)
    r,m=lsq(Z); print(f"  (r_sigma, r_pi) linear LSQ: RMSE {r:.3f} MaxAE {m:.2f}  CV-RMSE {lsq_cv(Z):.3f}")
    k=krr_cv(Z); print(f"  (r_sigma, r_pi) KRR CV-RMSE {k[0]:.3f} (lambda {k[1]}, gamma {k[2]})")
    a,l=sep(Z); print(f"  (r_sigma, r_pi) linear RS/ZB separation: {a*100:.1f}% fitted, {l*100:.1f}% leave-one-out")
    k=krr_cv(R4); print(f"  (rs_A, rp_A, rs_B, rp_B) KRR CV-RMSE {k[0]:.3f}")
    if name.startswith('this'):
        r,m=lsq(G); print(f"  Ghiringhelli 2D descriptor (with our IP, EA): LSQ RMSE {r:.3f} MaxAE {m:.2f}  CV-RMSE {lsq_cv(G):.3f}")
        r1,m1=lsq(G[:,:1]); print(f"  Ghiringhelli 1D descriptor: LSQ RMSE {r1:.3f} MaxAE {m1:.2f}")
        a,l=sep(G); print(f"  2D descriptor RS/ZB separation: {a*100:.1f}% fitted, {l*100:.1f}% LOO")
    print()
cmp=pd.concat([zun.add_suffix('_Zunger'),ours.add_suffix('_ours')],axis=1).loc[sorted(set(O.A)|set(O.B))]
print(cmp.round(3).to_string()); print('r(rs) %.3f  r(rp) %.3f'%(cmp.rs_Zunger.corr(cmp.rs_ours),cmp.rp_Zunger.corr(cmp.rp_ours)))
cmp.to_csv('radii_comparison_octet.csv')

import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Liberation Sans','font.size':28,'lines.linewidth':2,'axes.linewidth':2,'xtick.major.width':2,'ytick.major.width':2,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False,'legend.fontsize':22})
for tag,rad in (('zunger',zun),('ours',ours)):
    Z,_,_=feats(rad)
    fig,ax=plt.subplots(figsize=(8,8))
    ax.plot(Z[lab==0,0],Z[lab==0,1],'o',ms=10,mfc='none',mec='tab:red',mew=2,label='rock salt')
    ax.plot(Z[lab==1,0],Z[lab==1,1],'s',ms=10,color='tab:blue',label='zinc blende')
    Xs=(Z-Z.mean(0))/Z.std(0); m=LogisticRegression(C=1e4,max_iter=10000).fit(Xs,lab)
    xx=np.linspace(Z[:,0].min(),Z[:,0].max(),50); w=m.coef_[0]; b=m.intercept_[0]
    yy=((-b-w[0]*(xx-Z[:,0].mean())/Z[:,0].std())/w[1])*Z[:,1].std()+Z[:,1].mean()
    ax.plot(xx,yy,'k--')
    ax.set_xlabel(r'$r_\sigma$ (Å)'); ax.set_ylabel(r'$r_\pi$ (Å)')
    ax.set_ylim(Z[:,1].min()-0.05*np.ptp(Z[:,1]),Z[:,1].max()+0.25*np.ptp(Z[:,1])); ax.legend(loc='upper left',ncol=2,handletextpad=0.1,columnspacing=0.8)
    ax.set_title('Zunger 1980' if tag=='zunger' else 'this work',fontsize=26)
    fig.tight_layout(); fig.savefig(f'Fig_structmap_{tag}.png',dpi=130); fig.savefig(f'Fig_structmap_{tag}.eps')
