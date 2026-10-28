#!/usr/bin/env python3
"""Response-only audit of saved VAMPNet models; no training or MD changes."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/r2c3_mpl')
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
C2 = OUT.parent / 'R2C2_state_remodeling'
PACKAGE = C2 / 'response_package'
SAM = ROOT / 'revision_2026_09_30/state_sampling'
AUD = ROOT / 'revision_2026_09_30/state_audit'
sys.path.insert(0, str(AUD))
import audit_states as legacy
SPEC = importlib.util.spec_from_file_location('common_analysis', C2 / 'no_population_filter/analyze.py')
n = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(n)
s, p = n.s, n.p
INPUTS = set()
CHECKS = []
STRIDE = 10
SYSTEM_NAMES = {'FUS': 'FUS-LCD', 'A1': 'A1-LCD', 'TDP43': 'TDP-43-LCD'}


def read(path):
    path = Path(path)
    INPUTS.add(path)
    with path.open() as handle:
        return list(csv.DictReader(handle))


def write(name, rows):
    assert rows
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def protect():
    paths = set()
    for directory in ('project_figure', 'project_figure_set_Tl_Th',
                      'revision_2026_09_30/submitted_original/260822_JPCL/Submit'):
        paths.update(q for q in (ROOT / directory).rglob('*') if q.is_file()
                     and q.suffix in ('.png', '.pdf', '.svg', '.eps', '.py', '.csv', '.tex'))
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            for role in s.ROLES:
                d = s.model_dir(system, family, role)
                paths.update(d.glob('state_map_att*_raw_to_final.npy'))
                paths.update([d / 'avg_contacts.npy', d / 'avg_mu.npy'])
    return {str(q.relative_to(ROOT)): sha(q) for q in sorted(paths)}


def tag(system, family, role):
    return f'{system}_{family}_{s.TEMPS[system][role]}'


def base(system, family, role):
    return dict(system=system, family=family, role=role, temperature=s.TEMPS[system][role])


def hard_scores(a, b, nstates=5):
    joint = np.bincount(a * nstates + b, minlength=nstates**2).reshape(nstates, nstates)
    den = joint.sum(0) + joint.sum(1)
    dice = np.divide(2 * np.diag(joint), den, out=np.full(nstates, np.nan), where=den > 0)
    ii, jj = linear_sum_assignment(-joint)
    return joint, dice, joint[ii, jj].sum() / joint.sum()


def self_tests():
    a = np.array([0, 0, 1, 1])
    b = 1 - a
    joint, dice, optimum = hard_scores(a, b, 3)
    assert joint.trace() == 0 and optimum == 1
    assert np.isnan(dice[2])
    np.testing.assert_equal(dice[:2], 0)
    n.tests()


def condition(system, family, role):
    info = base(system, family, role)
    name = tag(system, family, role)
    stats_path = SAM / f'{name}_sufficient_statistics.npz'
    INPUTS.add(stats_path)
    with np.load(stats_path) as z:
        counts = z['counts'].sum((1, 2))
        masses = z['masses'].sum((1, 2))
        lengths, attempts = z['lengths'], z['attempts']
        assert int(z['lag']) == 200 and float(z['dt']) == .01
    # Corrected contact indices did not change any within-condition seed permutation.
    if system == 'TDP43':
        fixed_path = C2 / 'tdp43_index_fix' / name / 'fixed_alignment.npz'
        INPUTS.add(fixed_path)
        with np.load(fixed_path) as z:
            np.testing.assert_array_equal(z['attempts'], attempts)
            np.testing.assert_array_equal(z['new_to_old'], np.tile(np.arange(5), (20, 1)))
    model_t, model_mu = p.transition(counts), p.pm(masses)
    avg_t, avg_mu = model_t.mean(0), model_mu.mean(0)
    states = legacy.original.transition_states(system)
    old = legacy.load_condition(system, family, role, states)
    INPUTS.update(old['paths'])
    np.testing.assert_allclose(model_t, old['local_t'], atol=2e-14, rtol=0)
    np.testing.assert_allclose(masses / lengths.sum(), old['local_p'], atol=2e-14, rtol=0)
    # Saved float32 softmax rows sum to one up to rounding; the current estimator
    # explicitly normalizes masses, unlike the historical frame-count average.
    np.testing.assert_allclose(model_mu, old['local_p'], atol=1e-7, rtol=0)
    d = s.model_dir(system, family, role)
    inter = d / s.intermediate_name(system)
    qmodels, inventory = [], []
    for att in attempts:
        mp = d / f'state_map_att{int(att)}_raw_to_final.npy'
        INPUTS.add(mp)
        perm = np.argsort(np.load(mp).astype(int))
        qs = []
        qdir = inter / f'qcache_combined_n5_att{att}_lag50'
        for ti, length in enumerate(lengths):
            path = qdir / f'traj{ti}.npy'
            raw = np.load(path, mmap_mode='r')
            assert raw.shape == (length, 5)
            q = np.asarray(raw[::STRIDE][:, perm], dtype=np.float32)
            assert np.isfinite(q).all() and q.min() >= 0
            np.testing.assert_allclose(q.sum(-1), 1, atol=1e-5, rtol=0)
            qs.append(q)
            # Digest the exact sampled input, alongside file metadata, not multi-GB full caches.
            inventory.append(dict(**info, attempt=int(att), trajectory=ti, frames=int(length),
                                  sampled_frames=len(q), path=str(path), bytes=path.stat().st_size,
                                  mtime_ns=path.stat().st_mtime_ns,
                                  sampled_aligned_q_sha256=hashlib.sha256(q.tobytes()).hexdigest()))
        qmodels.append(np.concatenate(qs))
    q = np.stack(qmodels)
    hard = q.argmax(-1)
    sums = q.sum(1, dtype=np.float64)
    pairs, soft_values, hard_values = [], [], []
    for i in range(len(attempts)):
        for j in range(i + 1, len(attempts)):
            intersection = np.minimum(q[i], q[j]).sum(0, dtype=np.float64)
            soft = 2 * intersection / (sums[i] + sums[j])
            joint, dice, optimum = hard_scores(hard[i], hard[j])
            fixed_agreement = joint.trace() / joint.sum()
            assert optimum >= fixed_agreement - 1e-12
            soft_values.append(soft)
            hard_values.append(dice)
            pairs.append(dict(**info, attempt_a=int(attempts[i]), attempt_b=int(attempts[j]),
                              sampled_frames=q.shape[1], stride_frames=STRIDE,
                              hard_agreement=float(fixed_agreement),
                              optimized_hard_agreement=float(optimum),
                              soft_overlap=float(intersection.sum() / q.shape[1])))
    soft_values, hard_values = np.array(soft_values), np.array(hard_values)
    display = old['mapping']
    populations, agreement, models = [], [], []
    for k in range(5):
        h = hard_values[:, k]
        valid = h[np.isfinite(h)]
        populations.append(dict(**info, local_state=k+1,
            submitted_display_state=display[k]+1 if k in display else '',
            displayed_in_submitted_network=k in display, mean_population=float(avg_mu[k]),
            population_model_min=float(model_mu[:, k].min()),
            population_model_max=float(model_mu[:, k].max()),
            models_gt_1pct=int((model_mu[:, k] > .01).sum())))
        agreement.append(dict(**info, local_state=k+1, mean_population=float(avg_mu[k]),
            soft_dice_mean=float(soft_values[:, k].mean()),
            soft_dice_p05=float(np.quantile(soft_values[:, k], .05)),
            soft_dice_p95=float(np.quantile(soft_values[:, k], .95)),
            hard_dice_mean=float(valid.mean()) if len(valid) else '',
            hard_dice_valid_pairs=len(valid), total_model_pairs=len(pairs),
            hard_assignment_models=int(np.any(hard == k, axis=1).sum())))
        for i, att in enumerate(attempts):
            models.append(dict(**info, attempt=int(att), local_state=k+1,
                               population=float(model_mu[i, k]),
                               dominant_local_state=int(model_mu[i].argmax())+1))
    dominant = np.bincount(model_mu.argmax(1), minlength=5)
    summary = dict(**info, models=len(attempts), trajectories=len(lengths), frames=int(lengths.sum()),
        sampled_frames=q.shape[1], n_mean_states_gt_1pct=int((avg_mu > .01).sum()),
        mean_per_model_states_gt_1pct=float((model_mu > .01).sum(1).mean()),
        displayed_population=float(avg_mu[list(display)].sum()),
        omitted_population=float(avg_mu[[k for k in range(5) if k not in display]].sum()),
        max_mean_population=float(avg_mu.max()), mean_model_max_population=float(model_mu.max(1).mean()),
        dominant_model_counts=json.dumps(dominant.tolist()),
        hard_agreement=float(np.mean([r['hard_agreement'] for r in pairs])),
        optimized_hard_agreement=float(np.mean([r['optimized_hard_agreement'] for r in pairs])),
        soft_overlap=float(np.mean([r['soft_overlap'] for r in pairs])))
    old_summary = next(r for r in read(SAM / 'condition_agreement.csv')
                       if (r['system'], r['family'], r['role']) == (system, family, role))
    np.testing.assert_allclose(summary['hard_agreement'], float(old_summary['hard_overlap_mean']), atol=1e-12)
    np.testing.assert_allclose(summary['soft_overlap'], float(old_summary['soft_overlap_mean']), atol=1e-10)
    matrices = [dict(**info, source_local=i+1, target_local=j+1,
                    population=float(avg_mu[i]), transition_probability=float(avg_t[i, j]), lag_ns=2.)
                for i in range(5) for j in range(5)]
    CHECKS.append(dict(condition=name, cache_transition_max_error=float(abs(model_t-old['local_t']).max()),
                       raw_cache_population_max_error=float(abs(masses/lengths.sum()-old['local_p']).max()),
                       normalized_population_shift=float(abs(model_mu-old['local_p']).max()),
                       row_sum_max_error=float(abs(avg_t.sum(1)-1).max()),
                       all_seed_permutations_unchanged=True))
    print(json.dumps(summary), flush=True)
    return populations, agreement, models, summary, matrices, pairs, inventory


def distances():
    assignments = read(PACKAGE / 'Table_R2C2_1_state_assignments.csv')
    c2rows = read(PACKAGE / 'Table_R2C2_2_metrics.csv')
    c2contrasts = read(PACKAGE / 'Table_R2C2_3_contrasts.csv')
    cachepath = C2 / 'followup/bootstrap_draws.npz'
    INPUTS.add(cachepath)
    rows, contrasts, reproduction, mappings = [], [], [], []
    with np.load(cachepath) as z:
        for system in s.SYSTEMS:
            values = {}
            for family in s.FAMILIES:
                accepted = [r for r in assignments if r['system'] == system and r['family'] == family
                            and r['accepted'] == 'True']
                pairs = [(int(r['local_TL'])-1, int(r['local_TH'])-1) for r in accepted]
                for k, r in enumerate(accepted):
                    mappings.append(dict(common_pair=k+1, **r))
                for source in ('sampling', 'model_seeds', 'sampling_10ns', 'sampling_100ns'):
                    actual = 'sampling' if family == 'free' and source.startswith('sampling') else source
                    arrays = [(z[f'{system}_{family}_{role}_{actual}_t'],
                               z[f'{system}_{family}_{role}_{actual}_mu']) for role in s.ROLES]
                    d = n.independent_distance(arrays[0][0], arrays[1][0], arrays[0][1], arrays[1][1], pairs)
                    assert len(d) == 4001
                    values[family, source] = d
                    stat = p.summary(d)
                    if source in ('sampling', 'model_seeds'):
                        ref = next(r for r in c2rows if (r['system'],r['family'],r['uncertainty'],r['metric']) ==
                                   (system, family, source, 'D_common_pooled_weights'))
                        for key, v in stat.items():
                            np.testing.assert_allclose(v, float(ref[key]), atol=1e-14, rtol=0)
                    rows.append(dict(system=system, family=family, uncertainty=source,
                        n_common=len(pairs), lag_ns=2., bootstrap_replicates=4000, **stat))
                states = legacy.original.transition_states(system)
                old = [legacy.load_condition(system, family, role, states) for role in s.ROLES]
                tl, th = [x['t'].mean(0) for x in old]
                ml, mh = [x['p'].mean(0) for x in old]
                common = np.array([k in old[0]['retained'] & old[1]['retained'] for k in states])
                oldrows = read(AUD / 'distance_comparison.csv')
                for definition in ('submitted_edge_mask', 'unfiltered_display_union', 'common_submatrix'):
                    d = legacy.metric(definition, tl, th, ml, mh, common, old[0]['mask'], old[1]['mask'])
                    ref = next(r for r in oldrows if (r['system'],r['family'],r['definition']) ==
                               (system, family, definition))
                    np.testing.assert_allclose(d, float(ref['distance']), atol=1e-14, rtol=0)
                    reproduction.append(dict(system=system, family=family, definition=definition,
                        distance=float(d), status='historical audit only; not the revised estimator'))
                reproduction.append(dict(system=system, family=family,
                    definition='revised_common_all_states_no_population_prefilter',
                    distance=float(values[family,'model_seeds'][0]), status='revised full-cache estimate'))
            for source in ('sampling', 'model_seeds', 'sampling_10ns', 'sampling_100ns'):
                stat = p.summary(values['condensate', source] - values['free', source])
                if source in ('sampling', 'model_seeds'):
                    ref = next(r for r in c2contrasts if (r['system'],r['uncertainty'],r['metric']) ==
                               (system, source, 'D_common_pooled_weights'))
                    for key, value in stat.items():
                        np.testing.assert_allclose(value, float(ref[key]), atol=1e-14, rtol=0)
                contrasts.append(dict(system=system, uncertainty=source, **stat))
    write('Table_C3_4_temperature_distances.csv', rows)
    write('Table_C3_5_condensate_minus_free.csv', contrasts)
    write('Table_C3_6_number_reconciliation.csv', reproduction)
    write('common_state_correspondence.csv', mappings)
    wide = []
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            selected = {r['uncertainty']:r for r in rows if r['system']==system and r['family']==family}
            seed, sample = selected['model_seeds'], selected['sampling']
            wide.append(dict(system=system,family=family,n_common=seed['n_common'],
                full_cache_estimate=seed['estimate'],model_seed_ci_low=seed['ci_low'],model_seed_ci_high=seed['ci_high'],
                sampling_eligible_estimate=sample['estimate'],sampling_ci_low=sample['ci_low'],sampling_ci_high=sample['ci_high']))
    write('Table_C3_4_six_comparisons.csv',wide)
    return rows, contrasts


def plot(populations, agreements, summaries):
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':15, 'axes.labelsize':17,
                         'xtick.labelsize':16, 'ytick.labelsize':15, 'pdf.fonttype':42})
    pop = np.array([r['mean_population'] for r in populations], float).reshape(12,5)
    dice = np.array([r['soft_dice_mean'] for r in agreements], float).reshape(12,5)
    labels = [SYSTEM_NAMES[r['system']] + ' | ' + ('Free' if r['family']=='free' else 'Cond.') +
              ' | ' + (r'$T_{\mathrm{L}}$' if r['role']=='Tl' else r'$T_{\mathrm{H}}$') for r in summaries]
    fig, axes = plt.subplots(1,2,figsize=(13.5,9),sharey=True)
    for ax, array, lim, cmap, letter in zip(axes,[pop,dice],[.7,1.],['YlGnBu','cividis'],'AB'):
        im = ax.imshow(array, vmin=0, vmax=lim, cmap=cmap, aspect='auto')
        for i in range(12):
            for j in range(5):
                label = f'{100*array[i,j]:.1f}%' if letter=='A' else f'{array[i,j]:.2f}'
                rgb = im.cmap(im.norm(array[i,j]))[:3]
                lum = np.dot(rgb,[.2126,.7152,.0722])
                ax.text(j,i,label,ha='center',va='center',fontsize=14,color='black' if lum>.55 else 'white')
        for y in [3.5,7.5]:
            ax.axhline(y,color='white',lw=3)
        ax.set_xticks(range(5),[rf'$s_{{{j+1}}}$' for j in range(5)])
        ax.set_yticks(range(12),labels)
        ax.set_xlabel('Within-condition state')
        ax.text(-.06,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=24)
        cb = fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.06,pad=.12)
        cb.set_label(r'$\mu_i$' if letter=='A' else 'Soft assignment overlap')
        cb.set_ticks([0,.2,.4,.6] if letter=='A' else [0,.25,.5,.75,1])
    fig.subplots_adjust(left=.27,right=.96,top=.95,bottom=.10,wspace=.14)
    for ext in ('png','pdf'):
        fig.savefig(OUT / f'Fig_C3_state_population_and_agreement.{ext}',dpi=200)
    plt.close(fig)


def table_pdf(populations, matrices, distances_rows):
    with PdfPages(OUT/'Table_C3_complete_matrices.pdf') as pdf:
        for system in s.SYSTEMS:
            fig, axes = plt.subplots(2,2,figsize=(13,8.4))
            for ri,family in enumerate(s.FAMILIES):
                for ci,role in enumerate(s.ROLES):
                    ax = axes[ri,ci]
                    sel = [r for r in matrices if (r['system'],r['family'],r['role'])==(system,family,role)]
                    t = np.array([float(r['transition_probability']) for r in sel]).reshape(5,5)
                    mu = [float(r['mean_population']) for r in populations
                          if (r['system'],r['family'],r['role'])==(system,family,role)]
                    cells = [[f'{mu[k]:.5f}'] + [f'{v:.5f}' for v in t[k]] for k in range(5)]
                    ax.axis('off')
                    table = ax.table(cellText=cells,colLabels=[r'$\mu_i$']+[rf'$s_{{{k+1}}}$' for k in range(5)],
                                     rowLabels=[rf'$s_{{{k+1}}}$' for k in range(5)],
                                     cellLoc='center',bbox=[.03,.22,.97,.62])
                    table.auto_set_font_size(False)
                    table.set_fontsize(11)
                    for (row,col),cell in table.get_celld().items():
                        cell.set_linewidth(.5)
                        if row==0 or col==-1: cell.set_facecolor('#eeeeee')
                    env = 'Free chain' if family=='free' else 'Condensate'
                    temp = r'$T_{\mathrm{L}}$' if role=='Tl' else r'$T_{\mathrm{H}}$'
                    ax.set_title(f'{env} | {temp} | {s.TEMPS[system][role]}',fontsize=15)
                    r = next(r for r in distances_rows if (r['system'],r['family'],r['uncertainty'])==
                             (system,family,'model_seeds'))
                    ax.text(.5,.06,rf'Paired $D_T^{{(w)}}={r["estimate"]:.6f}$; '+f'{r["n_common"]} common states',
                            ha='center',transform=ax.transAxes,fontsize=12)
            fig.suptitle(SYSTEM_NAMES[system] + ' | Complete transition matrices at 2 ns',fontsize=19)
            fig.text(.5,.025,'Rows: origin; columns: destination. Local labels are condition-specific. '
                     'CSV files retain full precision.',ha='center',fontsize=10)
            fig.subplots_adjust(left=.05,right=.98,top=.89,bottom=.07,hspace=.19,wspace=.18)
            pdf.savefig(fig)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--render-only',action='store_true')
    args = parser.parse_args()
    self_tests()
    before = protect()
    if args.render_only:
        populations = read(OUT/'Table_C3_1_all_populations.csv')
        agreements = read(OUT/'Table_C3_2_state_agreement.csv')
        summaries = read(OUT/'condition_summary.csv')
        matrices = read(OUT/'Table_C3_3_full_transition_matrices.csv')
        drows = read(OUT/'Table_C3_4_temperature_distances.csv')
        for row in drows:
            row['estimate'] = float(row['estimate'])
    else:
        populations,agreements,models,summaries,matrices,pairs,inventory = [],[],[],[],[],[],[]
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                for role in s.ROLES:
                    out = condition(system,family,role)
                    for dest, data in zip([populations,agreements,models,matrices,pairs,inventory],
                                          [out[0],out[1],out[2],out[4],out[5],out[6]]):
                        dest.extend(data)
                    summaries.append(out[3])
        for filename, rows in [('Table_C3_1_all_populations.csv',populations),
                               ('Table_C3_2_state_agreement.csv',agreements),
                               ('Table_C3_3_full_transition_matrices.csv',matrices),
                               ('per_model_populations.csv',models),('condition_summary.csv',summaries),
                               ('all_2280_model_pairs.csv',pairs),('qcache_sample_inventory.csv',inventory)]:
            write(filename, rows)
        drows, contrasts = distances()
    plot(populations,agreements,summaries)
    table_pdf(populations,matrices,drows)
    after = protect()
    assert before == after, 'Protected original files changed during audit'
    if not args.render_only:
        (OUT/'protected_original_files.json').write_text(json.dumps(before,indent=2)+'\n')
        manifest = {str(path):sha(path) for path in sorted(INPUTS)}
        manifest[str(Path(__file__).resolve())] = sha(__file__)
        (OUT/'input_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
        (OUT/'validation.json').write_text(json.dumps(dict(conditions=CHECKS,
            protected_files_unchanged=len(before),sampled_q_files=len(inventory),
            model_pairs=len(pairs),states=len(populations),matrix_entries=len(matrices),
            matched_C2_estimates_and_intervals=True,tests_passed=True),indent=2)+'\n')
    print('Completed; original files unchanged:',len(before),flush=True)


if __name__ == '__main__':
    main()
