#!/usr/bin/env python3
"""Read-only R2C2 follow-up using saved VAMPNet sufficient statistics."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
SAM = ROOT / 'revision_2026_09_30/state_sampling'
AUD = ROOT / 'revision_2026_09_30/state_audit'
sys.path.insert(0, str(ROOT / 'project_figure/scripts/supplementary_figures'))
import shared as s
SPEC = importlib.util.spec_from_file_location('previous_results', OUT.parent / 'build_results.py')
old = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(old)
NBOOT = 4000
SEED = 20261008
INPUTS = set()
VALIDATION = []
FULL = old.read(AUD / 'full_local_transition_matrices.csv')
INPUTS.add(AUD / 'full_local_transition_matrices.csv')


def write(name, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def load_npz(path):
    INPUTS.add(path)
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def load_array(path):
    INPUTS.add(path)
    return np.load(path).astype(float)


def transition(c):
    c = .5 * (c + c.swapaxes(-1, -2))
    den = c.sum(-1, keepdims=True)
    assert np.all(den > 0), 'Unobserved source row; do not impute a zero or self row.'
    t = c / den
    np.testing.assert_allclose(t.sum(-1), 1., atol=1e-12)
    return t


def pm(m):
    assert np.all(m.sum(-1) > 0)
    return m / m.sum(-1, keepdims=True)


def draw_condition(system, family, role, rng):
    tag = f'{system}_{family}_{s.TEMPS[system][role]}'
    z = load_npz(SAM / f'{tag}_sufficient_statistics.npz')
    assert int(z['lag']) == 200 and np.isclose(z['dt'], .01)
    c, m = z['counts'].sum(2), z['masses'].sum(2)
    model_t, model_p = transition(c.sum(1)), pm(m.sum(1))
    reference = old.load_full(FULL, system, family, role)
    err_t = float(np.max(abs(model_t.mean(0) - reference['t'])))
    err_p = float(np.max(abs(model_p.mean(0) - reference['mu'])))
    assert err_t < 1e-8 and err_p < 1e-8, (tag, err_t, err_p)
    nm = len(model_t)
    w = rng.multinomial(nm, np.full(nm, 1 / nm), NBOOT) / nm
    result = {'model_seeds': (np.concatenate([model_t.mean(0)[None], np.einsum('bm,mij->bij', w, model_t)]),
                              np.concatenate([model_p.mean(0)[None], w @ model_p]))}
    excluded = []
    for block in ([0] if family == 'free' else [10, 50, 100]):
        if family == 'condensate':
            a = load_npz(SAM / f'{tag}_absolute_statistics.npz')
            factor = block // 10
            c = a['counts'].reshape(nm, -1, factor, 5, 5).sum(2)
            m = a['masses'].reshape(nm, -1, factor, 5).sum(2)
            excluded = a['excluded_trajectories'].tolist()
        nu = c.shape[1]
        w = rng.multinomial(nu, np.full(nu, 1 / nu), NBOOT)
        bt = transition(np.einsum('bu,muij->bmij', w, c, optimize=True)).mean(1)
        bp = pm(np.einsum('bu,mui->bmi', w, m, optimize=True)).mean(1)
        key = 'sampling' if block in (0, 50) else f'sampling_{block}ns'
        result[key] = (np.concatenate([transition(c.sum(1)).mean(0)[None], bt]),
                       np.concatenate([pm(m.sum(1)).mean(0)[None], bp]))
    VALIDATION.append(dict(condition=tag, models=nm, full_matrix_max_error=err_t,
                           full_population_max_error=err_p, excluded_trajectory_indices=excluded,
                           sampling_population_max_shift=float(abs(result['sampling'][1][0]-model_p.mean(0)).max())))
    return result


def summary(values):
    assert np.isfinite(values).all()
    lo, hi = np.quantile(values[1:], [.025, .975], axis=0)
    return dict(estimate=float(values[0]), ci_low=float(lo), ci_high=float(hi))


def maps(system, family):
    return [old.load_full(FULL, system, family, r)['mapping'] for r in s.ROLES]


def metrics(arrays, mapping):
    common = sorted(set(mapping[0]) & set(mapping[1]))
    ids = [[x[k] for k in common] for x in mapping]
    ts = [a[0][:, i][:, :, i] for a, i in zip(arrays, ids)]
    ps = [a[1][:, i] for a, i in zip(arrays, ids)]
    norm = [pm(p) for p in ps]
    w = .5 * (norm[0] + norm[1])
    delta = ts[1] - ts[0]
    out = {
        'population_TV_conditional_C': .5 * abs(norm[1] - norm[0]).sum(-1),
        'population_common_mass_change': ps[1].sum(-1) - ps[0].sum(-1),
        'D_common_pooled_weights': np.sqrt(np.sum(w * np.sum(delta**2, -1), -1)),
        'D_common_uniform_weights': np.sqrt(np.mean(np.sum(delta**2, -1), -1)),
        'D_common_fixed_TL_weights': np.sqrt(np.sum(norm[0] * np.sum(delta**2, -1), -1)),
        'D_common_fixed_TH_weights': np.sqrt(np.sum(norm[1] * np.sum(delta**2, -1), -1)),
    }
    mid = .5 * (norm[0] + norm[1])
    js = .5 * sum(np.sum(np.where(p > 0, p*np.log2(np.maximum(p, 1e-300)/mid), 0), -1) for p in norm)
    out['population_JS_distance_base2_conditional_C'] = np.sqrt(np.maximum(js, 0))
    return out, common, ids, ts, ps


def analyze_changes(data):
    rows, pop, entries, exits, contrast = [], [], [], [], []
    values = {}
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            mapping = maps(system, family)
            sources = data[system, family, 'Tl'].keys()
            for source in sources:
                arrays = [data[system, family, r][source] for r in s.ROLES]
                ms, common, ids, ts, ps = metrics(arrays, mapping)
                for metric, vals in ms.items():
                    values[system, family, source, metric] = vals
                    rows.append(dict(system=system, family=family, uncertainty=source,
                                     metric=metric, common_states=' '.join(str(k+1) for k in common), **summary(vals)))
                for a, state in enumerate(common):
                    base = dict(system=system, family=family, uncertainty=source, state=state+1)
                    pop.append(dict(**base, population_TL=float(ps[0][0,a]), population_TH=float(ps[1][0,a]),
                                    **summary(ps[1][:,a]-ps[0][:,a])))
                    exit_c = [1-t.sum(-1) for t in ts]
                    stay = [t[:,a,a] for t in ts]
                    for label, vals in [('exit_C', [x[:,a] for x in exit_c]),
                                        ('leave_state_at_2ns', [1-x for x in stay])]:
                        exits.append(dict(**base, quantity=label, value_TL=float(vals[0][0]), value_TH=float(vals[1][0]),
                                          **summary(vals[1]-vals[0])))
                    for b, target in enumerate(common):
                        entries.append(dict(**base, target_state=target+1, lag_ns=2.,
                                            probability_TL=float(ts[0][0,a,b]), probability_TH=float(ts[1][0,a,b]),
                                            **summary(ts[1][:,a,b]-ts[0][:,a,b])))
                # Probability conservation includes exit C, not a fabricated outside-source state.
                np.testing.assert_allclose((ts[1]-ts[0]).sum(-1)+(1-ts[1].sum(-1))-(1-ts[0].sum(-1)), 0, atol=1e-12)
            ref = old.one(old.read(OUT.parent / 'common_transition_distances_rechecked.csv'),
                          system=system, family=family, metric='common_submatrix')
            np.testing.assert_allclose(values[system,family,'model_seeds','D_common_pooled_weights'][0],
                                       float(ref['distance']), atol=1e-10)
    for system in s.SYSTEMS:
        for source in ['sampling', 'model_seeds']:
            for metric in ms:
                a = values[system, 'condensate', source, metric] - values[system, 'free', source, metric]
                contrast.append(dict(system=system, uncertainty=source, metric=metric, **summary(a)))
    write('population_changes.csv', pop)
    write('common_transition_entry_changes.csv', entries)
    write('state_exit_changes.csv', exits)
    write('separated_metrics.csv', rows)
    write('condensate_minus_free.csv', contrast)
    return pop, entries, exits, rows, contrast


def contact_features(a):
    i, j = np.triu_indices(a.shape[-1], 5)
    v = np.log10(np.clip(a[:, i, j], 1e-4, None))
    v -= v.mean(1, keepdims=True)
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def indexing_diagnostic(a):
    assert a.shape == (5, 146, 146)
    i, j = np.triu_indices(148, 3)
    x, y = np.triu_indices(146, 1)
    assert len(i) == len(x) == 10585
    v = a[:, x, y]
    result = np.zeros((5, 148, 148))
    result[:, i, j] = v
    result[:, j, i] = v
    np.testing.assert_array_equal(result[:, i, j], a[:, x, y])
    return result


def audit_matching(data):
    inventory, corr_rows, scans, diagnostic_distances = [], [], [], []
    corrected = {}
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            mapping = maps(system, family)
            common = set(mapping[0]) & set(mapping[1])
            contact = [load_array(s.model_dir(system, family, role) / 'avg_contacts.npy') for role in s.ROLES]
            for role, mapping_r, a in zip(s.ROLES, mapping, contact):
                d = s.model_dir(system, family, role)
                INPUTS.add(s.STATE_SOURCE[system][family])
                mu = load_array(d / 'avg_mu.npy')
                source_rows = old.read(s.STATE_SOURCE[system][family])
                for source_row in source_rows:
                    if source_row['temperature'] == s.TEMPS[system][role]:
                        assert abs(mu[int(source_row['local_state'])]-float(source_row['mu'])) < 1e-7
                inverse = {local: display for display, local in mapping_r.items()}
                for local in range(5):
                    display = inverse.get(local)
                    reason = ('retained_common' if display in common else
                              'retained_no_accepted_counterpart' if display is not None else 'omitted_by_population_selection')
                    inventory.append(dict(system=system, family=family, role=role, local_state=local,
                        display_state='' if display is None else display+1, population=float(mu[local]),
                        classification=reason, contact_nres=a.shape[-1],
                        selection_cutoff=.02 if system=='FUS' else .1,
                        cutoff_evidence='code_and_retained_set' if system!='TDP43' else 'consistent_with_retained_set_not_runtime_log',
                        interpretation='Omitted/unmatched does not demonstrate physical absence'))
            variants = {'saved_maps': contact}
            if system == 'TDP43':
                variants['index_reconstruction_diagnostic'] = [indexing_diagnostic(contact[0]), contact[1]]
                corrected[family] = variants['index_reconstruction_diagnostic'][0]
            for variant, (low, high) in variants.items():
                n = min(low.shape[-1], high.shape[-1])
                corr = contact_features(high[:,:n,:n]) @ contact_features(low[:,:n,:n]).T
                for hi in range(5):
                    for lo in range(5):
                        corr_rows.append(dict(system=system, family=family, variant=variant,
                                              local_TL=lo, local_TH=hi, log_contact_correlation=float(corr[hi,lo])))
                kept = [sorted(m.values()) for m in mapping]
                rr, cc = linear_sum_assignment(1-corr[np.ix_(kept[1], kept[0])])
                for cutoff in [.80,.85,.90,.95]:
                    accepted = [(kept[0][j], kept[1][i]) for i,j in zip(rr,cc)
                                if corr[kept[1][i],kept[0][j]] >= cutoff]
                    scans.append(dict(system=system, family=family, variant=variant, correlation_cutoff=cutoff,
                                      n_common=len(accepted), accepted_local_TL_TH=str(accepted),
                                      status='diagnostic only; FUS production uses a different matcher'))
                    if system == 'TDP43' and variant == 'index_reconstruction_diagnostic' and cutoff == .90:
                        newmap = [{k:l for k,(l,h) in enumerate(accepted)}, {k:h for k,(l,h) in enumerate(accepted)}]
                        for source in ['sampling','model_seeds']:
                            ms, *_ = metrics([data[system,family,r][source] for r in s.ROLES], newmap)
                            for metric, vals in ms.items():
                                diagnostic_distances.append(dict(system=system, family=family, uncertainty=source,
                                    metric=metric, common_pairs=str(accepted),
                                    status='not adopted; within-temperature seed alignment still fixed', **summary(vals)))
    write('state_inventory.csv', inventory)
    write('contact_match_correlations.csv', corr_rows)
    write('contact_match_threshold_diagnostic.csv', scans)
    write('TDP43_reindexed_mapping_sensitivity.csv', diagnostic_distances)
    np.savez_compressed(OUT / 'TDP43_contact_index_reconstruction_DIAGNOSTIC.npz', **corrected)
    return inventory


def save(fig, stem):
    for ext in ['png','pdf']:
        fig.savefig(OUT / f'{stem}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':13, 'axes.labelsize':14,
                         'axes.titlesize':15, 'xtick.labelsize':12, 'ytick.labelsize':12,
                         'pdf.fonttype':42, 'legend.frameon':False})


def plot_state_changes(pop, exits):
    style()
    fig, axs = plt.subplots(3,3,figsize=(15,10.5))
    for row, system in enumerate(s.SYSTEMS):
        for col, quantity in enumerate(['population', 'leave_state_at_2ns', 'exit_C']):
            ax = axs[row,col]
            table = pop if col==0 else [r for r in exits if r['quantity']==quantity]
            x, labels, k = [], [], 0
            for family in s.FAMILIES:
                subset = [r for r in table if r['system']==system and r['family']==family and r['uncertainty']=='sampling']
                for r in subset:
                    other = old.one(table, system=system, family=family, uncertainty='model_seeds', state=r['state'])
                    for offset, color, point in [(-.09,'#303030',r),(.09,'#397eb1',other)]:
                        ax.plot([k+offset]*2, 100*np.array([point['ci_low'],point['ci_high']]), color=color, lw=1.5)
                        ax.plot(k+offset, 100*point['estimate'], 'o', color=color, ms=5,
                                mfc='white' if family=='free' else color)
                    labels.append(r'$s_{%d}$'%r['state']);x.append(k);k+=1
                ax.text(np.mean(x[-len(subset):]), .98, 'Free chain' if family=='free' else 'Condensate',
                        transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=11)
                k+=.6
            ax.set_xticks(x, labels)
            ax.axhline(0, color='.6', ls=':', lw=1)
            ax.set_ylim([(-40,45),(-45,45),(-45,45)][col])
            if row==0:ax.set_title(['Population change','State-leaving change','Exit from common set'][col], pad=12)
            if col==0:ax.set_ylabel(s.DISPLAY[system]+'\nChange (percentage points)')
    fig.legend(handles=[Line2D([],[],color='#303030',marker='o',ls='none',label='Sampling'),
                        Line2D([],[],color='#397eb1',marker='o',ls='none',label='Model seeds')],
               loc='upper center',ncol=2,bbox_to_anchor=(.5,.96))
    fig.suptitle('Population and 2 ns transition changes',fontsize=20,y=.995)
    fig.subplots_adjust(left=.10,right=.985,top=.87,bottom=.10,wspace=.22,hspace=.26)
    fig.text(.5,.025,'High minus low temperature; pointwise 95% intervals. TDP-43 uses the original, provisional mapping.',
             ha='center',fontsize=11)
    save(fig,'C2_04_state_changes_PROVISIONAL')


def plot_transition_entries(entries):
    style()
    fig, axs = plt.subplots(3,2,figsize=(12,12))
    for i,system in enumerate(s.SYSTEMS):
        for j,family in enumerate(s.FAMILIES):
            ax=axs[i,j]
            subset=[r for r in entries if r['system']==system and r['family']==family and r['uncertainty']=='sampling']
            states=sorted({r['state'] for r in subset})
            a=np.array([[old.one(subset,state=x,target_state=y)['estimate'] for y in states] for x in states])
            im=ax.imshow(a,vmin=-.4,vmax=.4,cmap='RdBu_r',aspect='equal')
            for row,x in enumerate(states):
                for col,y in enumerate(states):
                    r=old.one(subset,state=x,target_state=y)
                    color='white' if abs(r['estimate'])>.25 else 'black'
                    ax.text(col,row,f"{r['estimate']:+.3f}\n[{r['ci_low']:+.3f},\n{r['ci_high']:+.3f}]",
                            ha='center',va='center',fontsize=10,color=color)
            ax.set_xticks(range(len(states)),[r'$s_{%d}$'%x for x in states])
            ax.set_yticks(range(len(states)),[r'$s_{%d}$'%x for x in states])
            ax.set_xlabel('To');ax.set_ylabel('From')
            ax.set_title(s.DISPLAY[system]+' | '+('Free chain' if family=='free' else 'Condensate'),pad=10)
    fig.subplots_adjust(left=.07,right=.89,bottom=.07,top=.91,hspace=.40,wspace=.22)
    cax=fig.add_axes([.92,.25,.018,.5]);fig.colorbar(im,cax=cax,label='Transition probability change')
    fig.suptitle('Common-state transition changes at 2 ns',fontsize=20,y=.98)
    fig.text(.5,.012,'High minus low; sampling 95% intervals. Original mapping; TDP-43 correspondence remains provisional.',
             ha='center',fontsize=10)
    save(fig,'C2_05_transition_changes_PROVISIONAL')


def plot_separated_metrics(rows):
    style()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.0))
    for ax, metric, title in zip(axes,
            ['population_JS_distance_base2_conditional_C', 'D_common_uniform_weights'],
            ['Population redistribution', 'Transition-matrix change']):
        for x, system in enumerate(s.SYSTEMS):
            for shift, family in [(-.14, 'free'), (.14, 'condensate')]:
                for offset, source, color in [(-.035, 'sampling', '#303030'), (.035, 'model_seeds', '#397eb1')]:
                    r = old.one(rows, system=system, family=family, uncertainty=source, metric=metric)
                    xx = x + shift + offset
                    ax.plot([xx]*2, [r['ci_low'], r['ci_high']], color=color, lw=1.5)
                    ax.plot(xx, r['estimate'], 'o', color=color, mfc='white' if family=='free' else color, ms=6)
        ax.set_xticks(range(3), ['FUS-LCD', 'A1-LCD', 'TDP-43-LCD'])
        ax.set_xlim(-.5,2.5);ax.set_ylim(0,.35);ax.set_title(title, pad=12)
        ax.set_ylabel('JS distance within C' if 'population' in metric else 'Uniform-weight distance within C')
    handles=[Line2D([],[],color='black',marker='o',mfc='white',ls='none',label='Free chain'),
             Line2D([],[],color='black',marker='o',ls='none',label='Condensate'),
             Line2D([],[],color='#303030',lw=2,label='Sampling'),
             Line2D([],[],color='#397eb1',lw=2,label='Model seeds')]
    fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(.5,.98), fontsize=11)
    fig.subplots_adjust(left=.09,right=.98,bottom=.18,top=.77,wspace=.30)
    fig.text(.5,.035,'Different statistics, not additive components. Original mapping; TDP-43 remains provisional.',
             ha='center',fontsize=10)
    save(fig,'C2_06_population_vs_transition_PROVISIONAL')


def main():
    rng=np.random.default_rng(SEED)
    data={}
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            for role in s.ROLES:
                print('Resampling',system,family,role,flush=True)
                data[system,family,role]=draw_condition(system,family,role,rng)
    pop, entries, exits, metrics_rows, contrasts=analyze_changes(data)
    audit_matching(data)
    plot_state_changes(pop,exits)
    plot_transition_entries(entries)
    plot_separated_metrics(metrics_rows)
    np.savez_compressed(OUT/'bootstrap_draws.npz', **{
        '_'.join((system,family,role,source,kind)):a
        for (system,family,role), bundle in data.items()
        for source, arrays in bundle.items() for kind,a in zip(['t','mu'],arrays)})
    INPUTS.update([Path(__file__), OUT.parent/'build_results.py', SAM/'analyze.py',
                   SAM/'absolute_time_exclusions.csv',
                   ROOT/'TDP43/tdp43_protein_single_simulation/thermo_figure3_state_population_contactmaps.py'])
    for family in s.FAMILIES:
        d = s.model_dir('TDP43', family, 'Tl')
        INPUTS.update([d/'logs/build_npz.log', d/'plot_fig3b.py', d/'scripts/build_npz.py'])
    manifest=[]
    for path in sorted(INPUTS):
        h=hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
        manifest.append(dict(path=str(path.relative_to(ROOT)),sha256=h.hexdigest(),bytes=path.stat().st_size))
    write('input_manifest.csv',manifest)
    (OUT/'validation.json').write_text(json.dumps(dict(bootstrap_replicates=NBOOT,seed=SEED,
        conditions=VALIDATION,lag_ns=2.,status='Original mapping conditional results; TDP43 indexing requires review',
        all_matrix_and_population_checks_passed=True),indent=2)+'\n')
    print('DONE:',OUT,flush=True)


if __name__=='__main__':
    main()
