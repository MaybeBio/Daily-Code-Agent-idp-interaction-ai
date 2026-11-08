import pickle,numpy as np,pandas as pd,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,sys
from idprg import distribution as D
SC=sys.argv[1]
real={r['name']:r for r in pickle.load(open('paper/cache/collected.pkl','rb'))['real']}
_z=np.load('runs/starling/benchmark/rg_samples.npz');R=pd.read_csv(SC+'/nu_check.csv').set_index('name');R['d']=(R.nu-R.nu_net)
A=R.reindex(R.d.abs().sort_values(ascending=False).index)
top=A.head(6);oth=A[(A.family!='A1_LCD')&(A.N>=20)&(~A.index.isin(top.index))].head(6)
fig,ax=plt.subplots(2,6,figsize=(15,5.2),sharex=True,sharey=True,layout='constrained')
x=np.linspace(.4,2.3,500);b=np.linspace(.4,2.3,46)
for row,(lab,T) in zip(ax,(('worst |Δν|, all',top),('worst |Δν|, no A1-LCD, N≥20',oth))):
    for a,(n,r) in zip(row,T.iterrows()):
        s=_z[n];p=real[n]['p'];s=s/s.mean();m=real[n]['mean']
        a.hist(s,bins=b,density=True,color='#2a78d6',alpha=.45,lw=0)
        a.plot(x,D.pdf(x*m,*p)*m,color='#D61431',lw=1.5)
        a.set_title(f"{n[:24]}  N={int(r.N)}",fontsize=8,loc='left')
        a.text(.97,.95,f"ν STARLING {r.nu:.2f}\nν network {r.nu_net:.2f}\nskew {r.skew_s:.2f} / {r.skew_n:.2f}\nW1 {r.w1csv:.2f}%",transform=a.transAxes,ha='right',va='top',fontsize=7.5)
        a.spines[['top','right']].set_visible(False)
    row[0].set_ylabel(lab+'\ndensity',fontsize=8)
for a in ax[1]:a.set_xlabel('Rg / <Rg>')
fig.savefig(SC+'/worst_nu.png',dpi=150)
print(top[['N','family','nu','nu_net','skew_s','skew_n','w1csv']].round(3).to_string());print(oth[['N','family','nu','nu_net','skew_s','skew_n','w1csv']].round(3).to_string())
print('rank of top by W1: median W1 of top 12 |dnu|',A.head(12).w1csv.median())
