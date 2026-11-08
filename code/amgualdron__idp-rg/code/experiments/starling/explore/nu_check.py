import pickle,numpy as np,pandas as pd,multiprocessing as mp
from scipy.stats import spearmanr
from idprg import mdfit,network,descriptors,distribution as D
from idprg.distribution import nu_of,u_of
real={r['name']:r for r in pickle.load(open('paper/cache/collected.pkl','rb'))['real']}
_z=np.load('runs/starling/benchmark/rg_samples.npz');S={k:_z[k] for k in _z.files};names=list(S)
t=pd.read_csv('runs/starling/benchmark/shape_w1.csv').set_index('name').loc[names]
net,spec=network.load();law=np.asarray(net.law,float)
d=descriptors.table(pd.Series([real[n] and pd.read_csv('data/benchmark.csv').set_index('name').sequence[n] for n in names]))
U=(np.arange(20000)+.5)/20000
def qshape(lnA,lnSm,nu):
    q=D.quantiles(U,lnA,lnSm,nu);return q/D.moments(lnA,lnSm,nu)[0]
def qs(x):x=np.sort(x);return x[(U*len(x)).astype(int)]/x.mean()
w=lambda a,b:100*np.abs(a-b).mean()
def one(i):
    n=names[i];r=real[n];x=np.asarray(S[n],float);xn=x/x.mean()
    lnA,lnSm,u,nll=mdfit.fit_free(xn);lnAr,_,ur,_=mdfit.fit_free(x)
    nu=float(nu_of(u));nun=float(r['p'][2]);lnAn=float(r['p'][0])
    uu,c=mdfit.compress(xn);N_=c.sum()
    L=lambda a,s,v:-np.sum(c*mdfit.logpdf(uu,a,s,v))/N_
    # nu fixed at network, A and scale refit
    a2,s2=mdfit._minimize(lambda t:L(t[0],t[1],nun),[[lnA,lnSm],[lnAn,lnSm]],[(-5,12),(None,None)])
    # A fixed at the network's, nu and scale refit
    s3,u3=mdfit._minimize(lambda t:L(lnAn,t[0],float(nu_of(t[1]))),[[lnSm,u],[lnSm,float(u_of(nun))]],[(None,None),(-8,8)])
    # A on the network's own law, nu and scale refit
    N=float(r['length']);h,sp=float(d.mean_hydro.iloc[i]),float(d.scd_pos.iloc[i])
    sl,ul=mdfit.fit_law(xn,N,h,sp,law,(lnSm,u))
    nul=float(nu_of(ul));Al=float(np.array([1,nul,nul**2,1/N,h,sp])@law)
    qx=qs(x);qf=qshape(lnA,lnSm,nu);qfix=qshape(a2,s2,nun);qn=qshape(lnAn,float(r['p'][1]),nun);qa=qshape(lnAn,s3,float(nu_of(u3)));ql=qshape(Al,sl,nul)
    return dict(name=n,nu=nu,nu_raw=float(nu_of(ur)),lnA=lnA,lnA_raw=lnAr,nu_net=nun,lnA_net=lnAn,a_fixnu=a2,nu_fixA=float(nu_of(u3)),nu_law=nul,lnA_law=Al,
       w_fit_samples=w(qf,qx),w_fix_samples=w(qfix,qx),w_net_samples=w(qn,qx),w_fit_fix=w(qf,qfix),w_fit_net=w(qf,qn),w_law_samples=w(ql,qx),w_fixA_samples=w(qa,qx),
       dll_fixnu=2*N_*(L(a2,s2,nun)-nll),dll_law=2*N_*(L(Al,sl,nul)-nll))
if __name__=='__main__':
    with mp.get_context('fork').Pool(12) as p: R=pd.DataFrame(p.map(one,range(len(names)))).set_index('name')
    R['N']=t.N;R['family']=t.family;R['skew_s']=t.skew_starling;R['skew_n']=t.skew_net;R['cv_s']=t.cv_starling;R['w1csv']=t.w1_net_starling
    R.to_csv('experiments/starling/explore/nu_check.csv')
    print('A normalise check: max|nu-nu_raw|',abs(R.nu-R.nu_raw).max(),'max|lnA diff|',abs(R.lnA-R.lnA_raw).max(), 'net check',abs(R.w_net_samples-R.w1csv).max())
    print('B fit vs samples W1 median',R.w_fit_samples.median(),' network vs samples',R.w_net_samples.median(),' floor .42')
    I=R[R.N>=20];A1=I[I.family=='A1_LCD'];O=I[I.family!='A1_LCD']
    print('C nu pinned at network, A refit: W1 to samples med',I.w_fix_samples.median(),' W1 to free fit med',I.w_fit_fix.median(),'LR med',I.dll_fixnu.median(),'frac LR>3.84',(I.dll_fixnu>3.84).mean(),'| A1 LR',A1.dll_fixnu.median(),'W1 fit-fix',A1.w_fit_fix.median(),'lnA fixnu med',I.a_fixnu.median())
    print('   law A: W1 to samples',I.w_law_samples.median(),'LR med',I.dll_law.median())
    for lab,d_ in (('all',I),('no A1',O),('A1',A1)):
        print(lab,len(d_),'| free nu rho %.2f mae %.3f | law-A nu rho %.2f mae %.3f | netA nu rho %.2f mae %.3f | lnA free med %.2f net %.2f rho %.2f'%(spearmanr(d_.nu,d_.nu_net).statistic,abs(d_.nu-d_.nu_net).mean(),spearmanr(d_.nu_law,d_.nu_net).statistic,abs(d_.nu_law-d_.nu_net).mean(),spearmanr(d_.nu_fixA,d_.nu_net).statistic,abs(d_.nu_fixA-d_.nu_net).mean(),d_.lnA.median(),d_.lnA_net.median(),spearmanr(d_.lnA,d_.lnA_net).statistic))
    print('nu vs lnA within free fits corr',np.corrcoef(I.nu,I.lnA)[0,1],'nu vs skew star rho',spearmanr(I.nu,I.skew_s).statistic,'nu vs cv rho',spearmanr(I.nu,I.cv_s).statistic)
    print('A1 skew star/net',A1.skew_s.median(),A1.skew_n.median(),' others',O.skew_s.median(),O.skew_n.median())
    print('W1(net,star) vs |dnu| rho',spearmanr(R.w1csv,abs(R.nu-R.nu_net)).statistic,' vs |dskew|',spearmanr(R.w1csv,abs(R.skew_s-R.skew_n)).statistic)
    # sensitivity: within the network's own predictions, shape W1 per difference in nu
    Q={n:qshape(*real[n]['p']) for n in names};P=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            P.append((abs(real[a]['p'][2]-real[b]['p'][2]),w(Q[a],Q[b])))
    P=np.array(P);
    for lo,hi in ((0,.01),(.01,.03),(.03,.06),(.06,.1),(.1,.2),(.2,.5)):
        k=(P[:,0]>=lo)&(P[:,0]<hi);print('net pairs |dnu| %.2f-%.2f n=%d median W1 %.2f'%(lo,hi,k.sum(),np.median(P[k,1])))
    # direct: free-fit shape at nu with the fitted A vs moved by 0.13 with A refit is w_fit_fix for A1
    print(A1[['nu','nu_net','nu_law','lnA','lnA_net','lnA_law','w_fit_samples','w_fit_net','w_fit_fix','dll_fixnu','skew_s','skew_n','w1csv']].round(3).head(6).to_string())
