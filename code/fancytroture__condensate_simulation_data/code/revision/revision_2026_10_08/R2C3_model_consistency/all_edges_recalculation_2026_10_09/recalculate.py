#!/usr/bin/env python3
"""Isolate edge-mask removal under the unchanged submitted display mapping."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/r2c3_all_edges_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT/'revision_2026_09_30/state_audit'))
import audit_states as old


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(name, data):
    with (OUT/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distance(tl, th, pl, ph):
    return np.sqrt(np.sum((pl+ph)[..., :, None]*.5*(th-tl)**2, axis=(-2,-1)))


def draws(data, indices):
    n = len(data['t'])
    # Seed-count weights preserve the original resampling while avoiding a
    # bootstrap-by-model-by-state-by-state intermediate array.
    weights = np.stack([(indices==i).sum(1) for i in range(n)], axis=1)/n
    return np.einsum('bm,mij->bij', weights, data['t']), weights @ data['p']


def evaluate(pair, indices=None, masked=False):
    if indices is None:
        values = [(d['t'].mean(0),d['p'].mean(0)) for d in pair]
    else:
        values = [draws(d, ids) for d,ids in zip(pair,indices)]
    tl,pl = values[0]
    th,ph = values[1]
    if masked:
        tl,th = tl*pair[0]['mask'],th*pair[1]['mask']
    return distance(tl,th,pl,ph)


def limits(values):
    return [float(x) for x in np.quantile(values,[.025,.975])]


def main():
    original = ROOT/'revision_2026_10_08/JCP_revision/submitted_original/260904_JCP_Communication'
    work = ROOT/'revision_2026_10_08/JCP_revision/working_copy'
    protected = [p for base in (original,work) for p in base.rglob('*') if p.is_file()]
    hashes = {str(p):digest(p) for p in protected}
    published = {(r['system'],r['family']):float(r['weighted_transition_distance']) for r in read(ROOT/'project_figure/Figure3_transition_distance_source.csv')}
    complete = read(OUT.parent/'Table_C3_3_full_transition_matrices.csv')
    common = read(OUT.parent/'Table_C3_4_six_comparisons.csv')
    results, contrasts, decomp, entries, state_mapping = [],[],[],[],[]
    checks, inputs, bootstrap = [],{},{}
    for path in (OUT.parent/'Table_C3_3_full_transition_matrices.csv',OUT.parent/'Table_C3_4_six_comparisons.csv',ROOT/'project_figure/Figure3_transition_distance_source.csv'):
        inputs[str(path)] = digest(path)
    point_seeds = {'FUS':3101,'A1':3201,'TDP43':3301}
    contrast_seeds = {'FUS':4101,'A1':4201,'TDP43':4301}
    for system in ('FUS','A1','TDP43'):
        states = old.original.transition_states(system)
        families = {}
        for family in ('free','condensate'):
            pair = [old.load_condition(system,family,role,states) for role in ('Tl','Th')]
            for path in (old.sh.STATE_SOURCE[system][family],old.original.RESULTS/f'transition_networks_{system.lower()}_free_condensate_global_aligned_source.csv'):
                inputs[str(path)] = digest(path)
            families[family] = pair
            for role,d in zip(('Tl','Th'),pair):
                assert len(d['t'])==20
                for path in d['paths']:
                    inputs[str(path)] = digest(path)
                local = d['local_t'].mean(0)
                supplied = np.zeros((5,5))
                for row in complete:
                    if (row['system'],row['family'],row['role'])==(system,family,role):
                        supplied[int(row['source_local'])-1,int(row['target_local'])-1] = float(row['transition_probability'])
                np.testing.assert_allclose(local,supplied,atol=2e-14,rtol=0)
                np.testing.assert_allclose(local.sum(1),1,atol=1e-12)
                for local_id,display_id in d['mapping'].items():
                    state_mapping.append(dict(system=system,family=family,role=role,local_state=local_id+1,display_state=display_id+1,population=float(d['local_p'][:,local_id].mean())))
                for i,si in enumerate(states):
                    for j,sj in enumerate(states):
                        exists = si in d['retained'] and sj in d['retained']
                        entries.append(dict(system=system,family=family,role=role,source_display=si+1,target_display=sj+1,measured_entry=exists,
                                            transition_probability=float(d['t'][:,i,j].mean()) if exists else '',
                                            legacy_zero_placeholder=not exists,old_edge_mask=int(d['mask'][i,j])))
            rng = np.random.default_rng(point_seeds[system]+(family=='condensate'))
            ids = rng.integers(0,20,(10000,2,20))
            d0 = float(evaluate(pair,masked=True))
            d1 = float(evaluate(pair))
            np.testing.assert_allclose(d0,published[system,family],atol=1e-7,rtol=0)
            for masked,label in [(True,'submitted_edge_mask'),(False,'all_entries_original_display_union')]:
                ds = evaluate(pair,[ids[:,0],ids[:,1]],masked=masked)
                lo,hi = limits(ds)
                point = d0 if masked else d1
                results.append(dict(system=system,family=family,definition=label,distance=point,model_ci_low=lo,model_ci_high=hi,n_bootstrap=10000))
                bootstrap[f'{system}_{family}_{label}'] = ds
            # Independently assemble the unmasked sum from local-state entries.
            local_sum = 0.
            reverse = [{v:k for k,v in d['mapping'].items()} for d in pair]
            means = [(d['local_t'].mean(0),d['local_p'].mean(0)) for d in pair]
            common_states = pair[0]['retained'] & pair[1]['retained']
            common_sum = 0.
            for source in states:
                p = [means[k][1][reverse[k][source]] if source in reverse[k] else 0. for k in (0,1)]
                for target in states:
                    ts = [means[k][0][reverse[k][source],reverse[k][target]] if source in reverse[k] and target in reverse[k] else 0. for k in (0,1)]
                    term = .5*sum(p)*(ts[1]-ts[0])**2
                    local_sum += term
                    if source in common_states and target in common_states:
                        common_sum += term
            np.testing.assert_allclose(np.sqrt(local_sum),d1,atol=2e-14)
            decomp.append(dict(system=system,family=family,squared_distance=local_sum,common_common_squared=common_sum,noncommon_squared=local_sum-common_sum,
                               noncommon_fraction=(local_sum-common_sum)/local_sum))
            c = next(r for r in common if (r['system'],r['family'])==(system,family))
            results.append(dict(system=system,family=family,definition='C2_common_no_population_prefilter',distance=float(c['full_cache_estimate']),model_ci_low=float(c['model_seed_ci_low']),model_ci_high=float(c['model_seed_ci_high']),n_bootstrap=4000))
            checks.append(dict(system=system,family=family,published_point_reproduced=True,full_matrix_export_matches=True,independent_local_sum_matches=True))
        rng = np.random.default_rng(contrast_seeds[system])
        ids = rng.integers(0,20,(20000,4,20))
        for masked,label in [(True,'submitted_edge_mask'),(False,'all_entries_original_display_union')]:
            free = evaluate(families['free'],[ids[:,0],ids[:,1]],masked)
            cond = evaluate(families['condensate'],[ids[:,2],ids[:,3]],masked)
            delta = cond-free
            lo,hi = limits(delta)
            contrasts.append(dict(system=system,definition=label,condensate_minus_free=float(evaluate(families['condensate'],masked=masked)-evaluate(families['free'],masked=masked)),
                                  model_ci_low=lo,model_ci_high=hi,n_bootstrap=20000))
            bootstrap[f'{system}_{label}_contrast'] = delta
    write('distance_comparison.csv',results)
    write('condensate_minus_free.csv',contrasts)
    write('squared_distance_decomposition.csv',decomp)
    write('original_display_mapping.csv',state_mapping)
    write('unmasked_original_display_matrices.csv',entries)
    write('complete_local_5_by_5_matrices.csv',complete)
    np.savez_compressed(OUT/'model_bootstrap_draws.npz',**bootstrap)
    assert hashes == {str(p):digest(p) for p in protected}
    assert all(digest(Path(p))==h for p,h in inputs.items())
    (OUT/'validation.json').write_text(json.dumps(dict(checks=checks,protected_files=len(hashes),all_protected_files_unchanged=True,
        all_inputs_unchanged=True,model_lag_ns=2,point_bootstrap=10000,contrast_bootstrap=20000,
        scope='Only edge-mask removal; original display support, alignment, populations, and zero padding unchanged; not a missing-state-free estimator'),indent=2)+'\n')
    (OUT/'input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
    render(results)
    print(json.dumps(contrasts,indent=2))


def render(results):
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'font.size':16,'axes.labelsize':17,'xtick.labelsize':16,'ytick.labelsize':16,'legend.fontsize':15.5})
    fig,axs = plt.subplots(1,2,figsize=(13,4.7),sharey=True)
    for ax,definition,title,letter in zip(axs,['submitted_edge_mask','all_entries_original_display_union'],['Original edge selection','No edge selection'],['A','B']):
        for family,offset,fill in [('free',-.12,'white'),('condensate',.12,'#4a4a4a')]:
            data = [next(r for r in results if (r['system'],r['family'],r['definition'])==(s,family,definition)) for s in ('FUS','A1','TDP43')]
            y = np.array([r['distance'] for r in data])
            ax.errorbar(np.arange(3)+offset,y,yerr=[y-[r['model_ci_low'] for r in data],np.array([r['model_ci_high'] for r in data])-y],fmt='o',ms=9,mfc=fill,mec='#4a4a4a',mew=1.5,color='#4a4a4a',capsize=4,lw=1.5)
        ax.set_xticks(range(3),['FUS-LCD','A1-LCD','TDP-43-LCD'])
        ax.set_ylim(0,.56)
        ax.set_yticks(np.arange(0,.6,.1))
        ax.set_xlim(-.5,2.5)
        ax.grid(axis='y',color='#ddd',alpha=.45)
        ax.set_title(title,fontsize=17,pad=15)
        ax.text(-.13,1.04,letter,transform=ax.transAxes,weight='bold',fontsize=22)
    axs[0].set_ylabel(r'$D_{\mathrm{T}}^{(w)}$')
    axs[0].legend(handles=[Line2D([],[],marker='o',linestyle='none',mfc='white',mec='#4a4a4a',ms=9,label='Free'),Line2D([],[],marker='o',linestyle='none',color='#4a4a4a',ms=9,label='Condensate')],loc='upper left',frameon=False)
    fig.subplots_adjust(left=.085,right=.985,bottom=.14,top=.86,wspace=.20)
    for suffix in ('png','pdf'):
        fig.savefig(OUT/f'Fig2B_edge_mask_only_comparison.{suffix}',dpi=220)
    plt.close(fig)


if __name__=='__main__':
    main()
