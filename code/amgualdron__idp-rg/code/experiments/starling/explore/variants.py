import sys,json,pathlib,numpy as np,pandas as pd,pickle,matplotlib
import make_figures as MF
from make_figures import plt,D,pdf,Line2D,MODEL,MOSS,DARK_ORANGE,INK,INK2,NET_STYLE,display_name,panel_label,COL2,spearmanr
SC=pathlib.Path(sys.argv[1]);MF.OUT=SC;MF.set_style()
real=pickle.load(open('paper/cache/collected.pkl','rb'))['real']
_z=np.load('runs/starling/benchmark/rg_samples.npz');S={k:np.asarray(_z[k],float) for k in _z.files}
t=pd.read_csv('runs/starling/benchmark/shape_w1.csv').set_index('name');fits=pd.read_csv('runs/starling/benchmark/nu_fits.csv').set_index('name')
rows={r['name']:r for r in real if r['name'] in S};names=list(rows);t=t.loc[names];t['nu_star']=fits.nu[names];t['nu_net']=[rows[n]['p'][2] for n in names]
def smp(x):
    x=x/x.mean();a,b,c=np.quantile(x,[.25,.5,.75]);return dict(bow=(c+a-2*b)/(c-a),pc=np.mean(x<.8),pe=np.mean(x>1.3),cv=x.std(ddof=1))
rng=np.random.default_rng(2026)
def blockboot(x,L=20,n=100):
    k=int(np.ceil(len(x)/L));out=[]
    for _ in range(n):
        st=rng.integers(0,len(x)-L+1,k);out.append(smp(np.concatenate([x[i:i+L] for i in st])[:len(x)]))
    return pd.DataFrame(out).std(ddof=1)
R=[]
for n in names:
    r=rows[n];p=r['p'];m=r['mean'];a,b,c=D.quantiles(np.array([.25,.5,.75]),*p)/m
    o=dict(name=n,bow_net=(c+a-2*b)/(c-a),pc_net=float(D.cdf(.8*m,*p)),pe_net=1-float(D.cdf(1.3*m,*p)),cv_net=t.cv_net[n])
    s,h=S[n],r['frames'];ss,hh=smp(s),smp(h)
    o|={k+'_star':v for k,v in ss.items()}|{k+'_hps1':v for k,v in hh.items()}
    e=pd.DataFrame([smp(rng.choice(s,len(s))) for _ in range(200)]).std(ddof=1);eh=blockboot(h)
    o|={k+'_se_star':e[k] for k in e.index}|{k+'_se_hps1':eh[k] for k in eh.index}
    R.append(o)
Q=pd.DataFrame(R).set_index('name');t=t.join(Q.drop(columns=[c for c in Q if c in t]))
I=t[t.N>=20];O=I[I.family!='A1_LCD'];A1=I[I.family=='A1_LCD']
M=[('cv','cv_net','cv_starling','cv_hps1','cv'),('bow','bow_net','bow_star','bow_hps1','bow'),('pc','pc_net','pc_star','pc_hps1','pc'),('pe','pe_net','pe_star','pe_hps1','pe')]
fams=[g for _,g in I.groupby('family')];boot=[pd.concat([fams[i] for i in rng.integers(len(fams),size=len(fams))]) for _ in range(2000)]
out=[];band={}
for lab,n_,s_,h_,e_ in M:
    w=1.96*I[e_+'_se_star'].median();band[lab]=w
    row=dict(metric=lab,sd_net=I[n_].std(),sd_star=I[s_].std(),sd_hps1=I[h_].std(),se_star=I[e_+'_se_star'].median(),se_hps1=I[e_+'_se_hps1'].median(),
      ceiling=np.mean([spearmanr(I[n_],I[n_]+rng.normal(0,I[e_+'_se_star'])).statistic for _ in range(200)]))
    for sn,d in (('all',I),('noA1',O),('A1',A1)):
        row[f'rho_star_{sn}']=spearmanr(d[n_],d[s_]).statistic;row[f'mae_star_{sn}']=(d[n_]-d[s_]).abs().mean();row[f'band_star_{sn}']=100*((d[n_]-d[s_]).abs()<=w).mean();row[f'bias_star_{sn}']=(d[n_]-d[s_]).mean()
    row['rho_hps1']=spearmanr(I[n_],I[h_]).statistic;row['mae_hps1']=(I[n_]-I[h_]).abs().mean();row['band_hps1']=100*((I[n_]-I[h_]).abs()<=w).mean()
    row['band_hps1_ownnoise']=100*((I[n_]-I[h_]).abs()<=1.96*I[e_+'_se_hps1'].median()).mean()
    row['rho_ci_lo'],row['rho_ci_hi']=np.percentile([spearmanr(g[n_],g[s_]).statistic for g in boot],[2.5,97.5])
    row['r2_vs_star']=1-((I[n_]-I[s_])**2).sum()/((I[s_]-I[s_].mean())**2).sum()
    out.append(row)
T=pd.DataFrame(out).set_index('metric').T;pd.set_option('display.width',200);print(T.round(4).to_string());T.to_csv(SC/'variants_table.csv')
print('rho bow vs nu star',spearmanr(I.bow_star,I.nu_star).statistic,'bow vs skew',spearmanr(I.bow_star,I.skew_starling).statistic,'pe vs nu',spearmanr(I.pe_star,I.nu_star).statistic,'pc vs nu',spearmanr(I.pc_star,I.nu_star).statistic,'pe vs skew',spearmanr(I.pe_star,I.skew_starling).statistic,'pc vs cv',spearmanr(I.pc_star,I.cv_starling).statistic,'pe vs cv',spearmanr(I.pe_star,I.cv_starling).statistic)
print('A1 medians bow star/net',A1.bow_star.median(),A1.bow_net.median(),'pc',A1.pc_star.median(),A1.pc_net.median(),'pe',A1.pe_star.median(),A1.pe_net.median());print('others bow',O.bow_star.median(),O.bow_net.median(),'pc',O.pc_star.median(),O.pc_net.median(),'pe',O.pe_star.median(),O.pe_net.median())
def scatter(ax,lab,name,xl,yl,ticks,spot,leg=False):
    n_,s_,h_=dict((m[0],m[1:4]) for m in M)[lab]
    vals=t[[n_,s_,h_]];lo,hi=vals.min().min(),vals.max().max();pad=.10*(hi-lo);lim=(lo-pad*(2.0 if spot[1]<.5 or leg else 1),hi+pad*(2.0 if spot[1]>.5 else 1))
    w=band[lab];fam=t.family=='A1_LCD';ins=t.N>=20
    ax.fill_between(lim,np.array(lim)-w,np.array(lim)+w,color='#d9d8d3',alpha=.75,linewidth=0,zorder=0);ax.plot(lim,lim,color=INK2,lw=.6,ls=(0,(4,3)),zorder=1)
    ax.scatter(t[h_][ins],t[n_][ins],s=8,marker='^',color=MOSS,alpha=.75,linewidths=0,zorder=2);ax.scatter(t[h_][~ins],t[n_][~ins],s=8,marker='^',facecolor='white',edgecolor=MOSS,linewidth=.6,zorder=2)
    ax.scatter(t[s_][ins&~fam],t[n_][ins&~fam],s=9,color=DARK_ORANGE,alpha=.7,linewidths=0,zorder=3);ax.scatter(t[s_][fam],t[n_][fam],s=9,color=MODEL,alpha=.8,linewidths=0,zorder=4);ax.scatter(t[s_][~ins],t[n_][~ins],s=9,facecolor='white',edgecolor=DARK_ORANGE,linewidth=.6,zorder=3)
    ax.set_xlim(lim);ax.set_ylim(lim);ax.set_aspect('equal',adjustable='box');ax.set_xticks(ticks);ax.set_yticks(ticks);ax.set_xlabel(xl);ax.set_ylabel(yl);ax.set_title(name,loc='left',pad=3,fontsize=7.5)
    c=T[lab];tx,ty,ha,va=spot
    ax.text(tx,ty,rf"STARLING $\rho$ {c.rho_star_all:.2f}, MAE {c.mae_star_all:.3f}"+"\n"+rf"HPS1 $\rho$ {c.rho_hps1:.2f}, MAE {c.mae_hps1:.3f}"+"\n"+rf"noise band $\pm${w:.3f}",transform=ax.transAxes,ha=ha,va=va,fontsize=5.6,linespacing=1.3,zorder=5)
    if leg:
        mk=lambda **k:Line2D([],[],ls='',ms=3,**k)
        ax.legend(handles=[mk(marker='o',color=MODEL,label='STARLING: A1-LCD'),mk(marker='o',color=DARK_ORANGE,label='STARLING: others'),mk(marker='^',color=MOSS,label='HPS1'),mk(marker='o',mfc='white',mec=INK2,label='N < 20')],loc='lower right',handletextpad=.1,borderaxespad=0,labelspacing=.1,fontsize=5.6)
def hists(axes,thr=False):
    o=t.w1_net_starling.sort_values();x=np.linspace(.4,2.2,400);b=np.linspace(.4,2.2,46);top=0
    for ax,pc,i in zip(axes,(25,50,75,100),range(4)):
        n=o.index[int(round(pc/100*(len(o)-1)))];r=rows[n];v=S[n]
        ax.hist(v/v.mean(),bins=b,density=True,color=MODEL,alpha=.45,histtype='stepfilled',linewidth=0,label='STARLING',zorder=1);ax.plot(x,pdf(x*r['mean'],r['p'])*r['mean'],label='network prediction',zorder=2,**NET_STYLE)
        if thr:
            for th in (.8,1.3): ax.axvline(th,color=INK2,lw=.6,ls=':',zorder=3)
        ax.set_xlim(.4,2.2);ax.set_title(display_name(n),loc='left',pad=3,fontsize=7);rank=f'largest of {len(o)}' if pc==100 else f'{pc}th percentile'
        ax.text(.97,.95,f"N = {int(t.N[n])}\n{rank}\n$W_1$ = {o[n]:.2f}%",transform=ax.transAxes,ha='right',va='top',fontsize=6,color=INK2,linespacing=1.3)
        ax.xaxis.set_major_locator(plt.MaxNLocator(4));ax.yaxis.set_major_locator(plt.MaxNLocator(4));ax.set_xlabel(r'$R_g/\langle R_g\rangle$');top=max(top,ax.get_ylim()[1])
    for ax in axes: ax.set_ylim(0,top*1.05)
    axes[0].set_ylabel('probability density');axes[0].legend(loc='center right',handlelength=1.3,borderaxespad=.2,fontsize=5.8)
# variant A
fig=plt.figure(figsize=(COL2,4.3),layout='constrained');outer=fig.add_gridspec(1,2,width_ratios=[1,2.05]);L=outer[0].subgridspec(2,1);Rr=outer[1].subgridspec(2,2)
a,b=fig.add_subplot(L[0]),fig.add_subplot(L[1]);h=[fig.add_subplot(Rr[i,j]) for i in (0,1) for j in (0,1)]
scatter(a,'cv','spread',r'reference $\sigma/\langle R_g\rangle$',r'network $\sigma/\langle R_g\rangle$',[.10,.15,.20,.25],(.97,.03,'right','bottom'))
scatter(b,'bow','quartile asymmetry (Bowley)','reference Bowley skewness','network Bowley skewness',[0,.1,.2],(.03,.97,'left','top'),leg=True)
hists(h)
for ax,l in zip([a,b]+h,'abcdef'): panel_label(ax,l,x=-.06 if l in 'cdef' else -.16,y=1.03)
MF.save(fig,'variant_A_bowley',tight=False)
# variant B
fig=plt.figure(figsize=(COL2,4.6),layout='constrained');outer=fig.add_gridspec(2,1,height_ratios=[1.25,1]);U=outer[0].subgridspec(1,3);Dn=outer[1].subgridspec(1,4)
a,b,c=[fig.add_subplot(U[0,i]) for i in range(3)];h=[fig.add_subplot(Dn[0,i],sharey=None) for i in range(4)]
scatter(a,'cv','spread',r'reference $\sigma/\langle R_g\rangle$',r'network $\sigma/\langle R_g\rangle$',[.10,.15,.20,.25],(.97,.03,'right','bottom'))
scatter(b,'pc',r'compact: $P(R_g<0.8\langle R_g\rangle)$','reference','network',[0,.1,.2,.3],(.97,.03,'right','bottom'))
scatter(c,'pe',r'extended: $P(R_g>1.3\langle R_g\rangle)$','reference','network',[0,.05,.10,.15],(.03,.97,'left','top'),leg=True)
hists(h,thr=True)
for ax,l in zip([a,b,c]+h,'abcdefg'): panel_label(ax,l,x=-.12,y=1.03)
MF.save(fig,'variant_B_tails',tight=False)
t[['N','family','bow_net','bow_star','bow_hps1','pc_net','pc_star','pc_hps1','pe_net','pe_star','pe_hps1','nu_star','nu_net']].to_csv(SC/'variants_values.csv')
