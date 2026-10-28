#!/usr/bin/env python3
"""Prepare response-only population/transition analysis on fixed all-state matches."""
import csv
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/r2c2_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('no_filter', OUT.parent/'no_population_filter/analyze.py')
n = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(n)
a, p, s = n.a, n.p, n.s
METRICS = ('population_TV_conditional_C', 'D_common_pooled_weights', 'D_common_uniform_weights',
           'D_common_fixed_TL_weights', 'D_common_fixed_TH_weights')
NAMES = ['FUS-LCD', 'A1-LCD', 'TDP-43-LCD']


def write(name, rows):
    assert rows
    with (OUT/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def stat(**fields):
    return fields


def validate_metrics():
    # Pure occupancy change at fixed transitions must not create a transition distance.
    t = np.array([[[.8,.2],[.3,.7]]])
    ml, mh = np.array([[.8,.2]]), np.array([[.3,.7]])
    mapping = [{0:0,1:1}, {0:0,1:1}]
    values = p.metrics([(t,ml),(t,mh)],mapping)[0]
    np.testing.assert_allclose(values['population_TV_conditional_C'], .5)
    for name in METRICS[1:]:
        np.testing.assert_allclose(values[name], 0.)
    # Fixed occupancy and a changed row isolate the transition calculation.
    th = np.array([[[.6,.4],[.3,.7]]])
    values = p.metrics([(t,ml),(th,ml)],mapping)[0]
    np.testing.assert_allclose(values['population_TV_conditional_C'], 0.)
    np.testing.assert_allclose(values['D_common_uniform_weights'], .2)
    np.testing.assert_allclose(values['D_common_pooled_weights'], np.sqrt(.064))
    # Unlike the distance above, zero-imputed missing rows yield a spurious term.
    assert np.linalg.norm(t[0,0]-np.zeros(2)) > 0
    n.tests()


def main_data():
    validate_metrics()
    nprotected = a.guard('after')
    distance, contrast, coverage, local_pops, common_pops, full_t, common_t = [], [], [], [], [], [], []
    assignments, draws, plot_data = [], {}, {}
    old_distance = a.read(OUT.parent/'no_population_filter/distances.csv')
    with np.load(OUT.parent/'followup/bootstrap_draws.npz') as cache:
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                cm, mu = zip(*(a.contacts(system,family,role) for role in s.ROLES))
                corr = p.contact_features(cm[1]) @ p.contact_features(cm[0]).T
                all_pairs, pairs = n.assign(corr)
                plot_data[system,family] = (corr,mu,pairs)
                maps = [{i:x[ri] for i,x in enumerate(pairs)} for ri in (0,1)]
                for l,h in all_pairs:
                    assignments.append(stat(system=system,family=family,local_TL=l+1,local_TH=h+1,
                        population_TL=float(mu[0][l]),population_TH=float(mu[1][h]),
                        contact_correlation=float(corr[h,l]),accepted=(l,h) in pairs))
                for source in ('sampling','model_seeds'):
                    arrays = [(cache[f'{system}_{family}_{role}_{source}_t'],
                               cache[f'{system}_{family}_{role}_{source}_mu']) for role in s.ROLES]
                    vals,_,ids,ts,ps = p.metrics(arrays,maps)
                    norms = [p.pm(x) for x in ps]
                    independent_tv = .5*np.abs(norms[1]-norms[0]).sum(-1)
                    independent_unif = np.linalg.norm(ts[1]-ts[0],axis=(1,2))/np.sqrt(len(pairs))
                    np.testing.assert_allclose(vals[METRICS[0]], independent_tv,atol=1e-14)
                    np.testing.assert_allclose(vals[METRICS[2]], independent_unif,atol=1e-14)
                    old = next(r for r in old_distance if (r['system'],r['family'],r['uncertainty'])==
                               (system,family,source))
                    for key,value in p.summary(vals[METRICS[1]]).items():
                        np.testing.assert_allclose(value,float(old[key]),atol=1e-14)
                    for metric in METRICS:
                        draws[system,family,source,metric] = vals[metric]
                        distance.append(stat(system=system,family=family,uncertainty=source,metric=metric,
                                             n_common=len(pairs),**p.summary(vals[metric])))
                    for ri,role in enumerate(s.ROLES):
                        t,m = arrays[ri]
                        np.testing.assert_allclose(t.sum(-1),1.,atol=1e-10)
                        np.testing.assert_allclose(m.sum(-1),1.,atol=1e-10)
                        mass = ps[ri].sum(-1)
                        other = m[:,[i for i in range(5) if i not in ids[ri]]].sum(-1)
                        np.testing.assert_allclose(mass+other,1.,atol=1e-10)
                        coverage.append(stat(system=system,family=family,role=role,uncertainty=source,
                                             **p.summary(mass)))
                        for local in range(5):
                            local_pops.append(stat(system=system,family=family,role=role,uncertainty=source,
                                local_state=local+1,common=local in ids[ri],**p.summary(m[:,local])))
                        if source=='model_seeds':
                            for i in range(5):
                                for j in range(5):
                                    full_t.append(stat(system=system,family=family,role=role,source_local=i+1,
                                                       target_local=j+1,T=float(t[0,i,j])))
                    for i,(l,h) in enumerate(pairs):
                        for kind,delta in [('raw_population',ps[1][:,i]-ps[0][:,i]),
                                           ('conditional_common_population',norms[1][:,i]-norms[0][:,i])]:
                            common_pops.append(stat(system=system,family=family,uncertainty=source,
                                local_TL=l+1,local_TH=h+1,kind=kind,**p.summary(delta)))
                        for j,(lj,hj) in enumerate(pairs):
                            common_t.append(stat(system=system,family=family,uncertainty=source,
                                source_TL=l+1,source_TH=h+1,target_TL=lj+1,target_TH=hj+1,
                                **p.summary(ts[1][:,i,j]-ts[0][:,i,j])))
                    np.testing.assert_allclose((norms[1]-norms[0]).sum(-1),0.,atol=1e-12)
                    assert np.all(vals[METRICS[0]]>=0) and np.all(vals[METRICS[0]]<=1)
            for source in ('sampling','model_seeds'):
                for metric in METRICS:
                    delta = draws[system,'condensate',source,metric]-draws[system,'free',source,metric]
                    contrast.append(stat(system=system,uncertainty=source,metric=metric,**p.summary(delta)))
    for filename,rows in [('Table_R2C2_1_state_assignments.csv',assignments),
                          ('Table_R2C2_2_metrics.csv',distance),('Table_R2C2_3_contrasts.csv',contrast),
                          ('Table_R2C2_4_common_coverage.csv',coverage),('all_local_populations.csv',local_pops),
                          ('common_population_changes.csv',common_pops),('full_transition_matrices.csv',full_t),
                          ('common_transition_changes.csv',common_t)]:
        write(filename,rows)
    np.savez_compressed(OUT/'metric_bootstrap_draws.npz',**{'__'.join(k):v for k,v in draws.items()})
    return distance,contrast,coverage,plot_data,nprotected


def interval(ax,x,r,marker='o',fill='white',size=8,color='.2',edgewidth=1.):
    lo,hi = r['ci_low'],r['ci_high']
    ax.vlines(x,lo,hi,color=color,lw=1.5)
    ax.hlines([lo,hi],x-.035,x+.035,color=color,lw=1.5)
    ax.plot(x,r['estimate'],marker=marker,mfc=fill,mec=color,ms=size,mew=edgewidth)


def legend(fig,y):
    fig.legend(handles=[Line2D([],[],marker='o',mfc='white',mec='.2',ls='',ms=8,label='Free chain'),
                        Line2D([],[],marker='s',color='.2',ls='',ms=8,label='Condensate')],
               loc='upper center',bbox_to_anchor=(.5,y),ncol=2,frameon=False)


def save(fig,name):
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'{name}.{ext}',dpi=220)
    plt.close(fig)


def metric_plot(distance,contrast,metrics,labels,name,title,with_contrast=False,
                sources=('sampling','model_seeds')):
    single_row = len(sources)==1
    fig,axes = plt.subplots(len(sources),2,figsize=(11.4,4.9 if single_row else 8),
                            sharey='col',squeeze=False)
    for ri,source in enumerate(sources):
        for ci,ax in enumerate(axes[ri]):
            metric = metrics[0] if with_contrast else metrics[ci]
            for i,system in enumerate(s.SYSTEMS):
                if with_contrast and ci==1:
                    r = next(r for r in contrast if (r['system'],r['uncertainty'],r['metric'])==
                             (system,source,metric))
                    interval(ax,i,r,fill='.2')
                else:
                    for fam,dx,mark,fill in [('free',-.1,'o','white'),('condensate',.1,'s','.2')]:
                        r = next(r for r in distance if (r['system'],r['family'],r['uncertainty'],r['metric'])==
                                 (system,fam,source,metric))
                        interval(ax,i+dx,r,mark,fill)
            ax.set_xticks(range(3),NAMES)
            ax.set_xlim(-.45,2.45)
            ax.set_ylabel(labels[ci])
            ax.set_title('Trajectory sampling' if source=='sampling' else 'Model seeds')
            ax.text(-.12,1.045,'ABCD'[ri*2+ci],transform=ax.transAxes,fontweight='bold',fontsize=15)
            ax.grid(axis='y',alpha=.18)
            ax.set_axisbelow(True)
            if with_contrast and ci==1:
                ax.axhline(0,color='.5',ls=':',lw=1.2)
    for ci in range(2):
        metric = metrics[0] if with_contrast else metrics[ci]
        rows = [r for r in (contrast if with_contrast and ci==1 else distance)
                if r['metric']==metric and r['uncertainty'] in sources]
        lo = min(r['ci_low'] for r in rows)
        hi = max(r['ci_high'] for r in rows)
        pad = .08*(hi-min(0,lo))
        axes[0,ci].set_ylim(min(0,lo)-pad if with_contrast and ci==1 else 0,hi+pad)
    fig.suptitle(title,fontsize=18,y=.98)
    legend(fig,.91 if single_row else .94)
    fig.subplots_adjust(top=.73 if single_row else .82,bottom=.14 if single_row else .08,
                        left=.09,right=.98,hspace=.40,wspace=.36)
    save(fig,name)


def matching_plot(plot_data):
    fig,axes = plt.subplots(3,2,figsize=(13,14))
    for ri,system in enumerate(s.SYSTEMS):
        for ci,family in enumerate(s.FAMILIES):
            ax = axes[ri,ci]
            corr,mu,pairs = plot_data[system,family]
            im = ax.imshow(corr,cmap='cividis',vmin=0,vmax=1,aspect='equal')
            for h in range(5):
                for l in range(5):
                    ax.text(l,h,f'{corr[h,l]:.2f}',ha='center',va='center',fontsize=11,
                            color='black' if corr[h,l]>.58 else 'white')
            for l,h in pairs:
                ax.add_patch(Rectangle((l-.46,h-.46),.92,.92,fill=False,ec='black',lw=2))
            ax.set_xticks(range(5),[rf'$s_{{{k+1}}}$'+'\n'+f'{100*mu[0][k]:.1f}%' for k in range(5)])
            ax.set_yticks(range(5),[rf'$s_{{{k+1}}}$'+f'  {100*mu[1][k]:.2f}%' for k in range(5)])
            ax.set_xlabel(r'$T_{\mathrm{L}}$ local state')
            ax.set_ylabel(r'$T_{\mathrm{H}}$ local state')
            family_label = 'Free chain' if family=='free' else 'Condensate'
            ax.set_title(f'{NAMES[ri]} | {family_label}',fontsize=15,pad=10)
            ax.text(-.32,1.10,'ABCDEF'[ri*2+ci],transform=ax.transAxes,fontweight='bold',fontsize=15)
            ax.tick_params(labelsize=11)
    fig.suptitle('State correspondence across temperatures',fontsize=19,y=.985)
    fig.subplots_adjust(left=.12,right=.97,bottom=.115,top=.915,hspace=.53,wspace=.50)
    cax = fig.add_axes([.25,.040,.52,.012])
    fig.colorbar(im,cax=cax,orientation='horizontal',label='Contact-map correlation')
    save(fig,'Fig_R2C2_1_state_correspondence')


def combined_plot(distance):
    fig,axes = plt.subplots(1,2,figsize=(11.4,4.3))
    sizes = {'free':np.sqrt(95),'condensate':np.sqrt(105)}
    color = '#4a4a4a'
    for ci,(ax,metric,label) in enumerate(zip(axes,[METRICS[1],METRICS[0]],
                                           [r'$D_T^{(w)}$',r'$D_{\mu,C}$'])):
        rows = [r for r in distance if r['metric']==metric and r['uncertainty']=='sampling']
        for i,system in enumerate(s.SYSTEMS):
            for family,dx,fill in [('free',-.1,'white'),('condensate',.1,color)]:
                row = next(r for r in rows if r['system']==system and r['family']==family)
                interval(ax,i+dx,row,'o',fill,size=sizes[family],color=color,edgewidth=1.8)
        ax.set_xticks(range(3),NAMES)
        ax.set_xlim(-.45,2.45)
        ax.set_ylim(0,1.08*max(r['ci_high'] for r in rows))
        ax.set_ylabel(label,fontsize=17,labelpad=10)
        ax.tick_params(axis='both',labelsize=16)
        ax.text(-.12,1.02,'AB'[ci],transform=ax.transAxes,fontweight='bold',fontsize=24,va='bottom')
        ax.grid(axis='y',alpha=.18)
        ax.set_axisbelow(True)
    axes[0].set_ylim(axes[1].get_ylim())
    for ax in axes:
        ax.set_yticks(np.arange(0.,.41,.1))
    axes[0].legend(handles=[
        Line2D([],[],marker='o',ls='',ms=sizes['free'],mfc='white',mec=color,mew=1.8,label='Free'),
        Line2D([],[],marker='o',ls='',ms=sizes['condensate'],mfc=color,mec=color,mew=1.8,label='Condensate')],
        loc='upper left',ncol=2,fontsize=15.5,frameon=False,
        handletextpad=.45,columnspacing=1.,borderaxespad=.3)
    fig.subplots_adjust(top=.88,bottom=.15,left=.09,right=.98,wspace=.36)
    save(fig,'Fig_R2C2_common_transition_and_population')


def main():
    plt.rcParams.update({'font.size':13,'axes.labelsize':14,'axes.titlesize':15,
                         'xtick.labelsize':12,'ytick.labelsize':12,'pdf.fonttype':42,'ps.fonttype':42})
    distance,contrast,coverage,plot_data,nprotected = main_data()
    combined_plot(distance)
    assert nprotected == a.guard('after')
    base = json.loads((OUT.parent/'no_population_filter/validation_and_settings.json').read_text())
    base.update(package_script_sha256=a.sha(Path(__file__)),
        population_distance='0.5*sum_C(abs(mu_H/sum_C(mu_H)-mu_L/sum_C(mu_L)))',
        uniform_transition_distance='sqrt(sum_C,sum_C(delta_T**2)/K), K=number of matched pairs',
        additional_tests=['Occupancy-only change gives zero transition distances',
                          'Transition-only change gives zero population TV',
                          'Independent TV/Frobenius and previous weighted-distance reproduction for all draws',
                          'All 12 full transition matrices stochastic; populations sum to 1',
                          'Common plus noncommon mass sums to 1; conditional population changes sum to 0',
                          f'{nprotected} original-file SHA256 checks passed'],
        original_files_unchanged=True,files_sha256={x.name:a.sha(x) for x in OUT.iterdir()
                                                 if x.suffix in ('.csv','.png','.pdf','.npz')})
    (OUT/'validation.json').write_text(json.dumps(base,indent=2)+'\n')
    print('Analysis, combined two-panel figure and validation complete.',flush=True)
    for row in contrast:
        if row['uncertainty']=='model_seeds': print(row,flush=True)


if __name__=='__main__':
    main()
