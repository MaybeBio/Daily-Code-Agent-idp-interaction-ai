"""
phospho_nmr_page4_prediction.py — the model's prediction for PAGE4 T51 phosphorylation from BENDER-MT's EXISTING
CALVADOS-2 runs (page4_CALVADOS, WT and HIPK1T51, 150 mM, pH 7.5, one run each, 278 K and 300 K). Written, hashed and run
before any chemical shift difference of BMRB 26544/26545 was computed (phospho_nmr_plan.md, hash 1851cf7c: this is the
plan's PAGE4 secondary, computed early so that the prediction is fixed before the data are seen). Descriptive.

Per temperature, on 1,000 evenly spaced frames of the second half (bender_common.load_calvados, as BENDER-MT):
  Rg change; contact frequencies p_ij (Calpha < 0.8 nm, |i - j| >= 3); per residue the net contact change
  sum_j dp_ij and the rewiring r_i = sum_j |dp_ij| (the plan's single-pair r_i); a rough noise floor from the two halves of
  the 1,000 frames of each run (within-run halves are autocorrelated, so the floor is a lower bound); the contacts of pT51
  and of the basic N-motif (residues 4-12) with the rest of the chain; the largest single contact changes.
At 300 K also the charge-model reweighting of the WT run (dq = -2 at T51, phospho_network.page4_phi, BETA at 300 K),
compared with the direct profile.

  BENDER_REANALYSIS_DIR=$BSS_REPO/bender_reanalysis ~/.venvs/finches/bin/python phospho_nmr_page4_prediction.py
      ->  phospho_nmr/page4_prediction_existing_runs.{txt,npz,png}
"""
import os
import sys

os.environ.setdefault('BENDER_REANALYSIS_DIR', (__import__('os').environ['BSS_REPO'] + '/bender_reanalysis'))
DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)

import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import phospho_network as pn  # noqa: E402

BASE = (__import__('os').environ['BSS_DATA'] + '/BENDER-MT/page4_CALVADOS')
OUT = os.path.join(DYN, 'phospho_nmr')
SITE = 51
NMOTIF = range(4, 13)
SEQ = 'MSARVRSRSRGRGDGQEAPDVVAFVAPGESQQEEPPTDNQDIEPGQEREGTPPIEERKVEGDCQEMDLEKTRSERGDGSDVKEKTPPNPKHAKTKEAGDGQP'


def load(T, form):
    d = os.path.join(BASE, T, f'PAGE4-{form}-I150-pH75')
    xyz, seq, info = pn.sp.bc.load_calvados(os.path.join(d, f'PAGE4-{form}-I150-pH75.dcd'), os.path.join(d, 'top.pdb'),
                                            k=1000)
    return np.asarray(xyz, float), info


def contacts(xyz):
    n = xyz.shape[1]
    far = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    C = np.zeros((len(xyz), n, n), bool)
    for f in range(len(xyz)):
        D = np.linalg.norm(xyz[f][:, None] - xyz[f][None], axis=2)
        C[f] = (D < 0.8) & far
    return C


def rg(xyz):
    return np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))


def main():
    os.makedirs(OUT, exist_ok=True)
    L = ['PAGE4 T51 phosphorylation: prediction from BENDER-MT existing CALVADOS-2 runs (one run per form; 150 mM, pH 7.5)',
         'computed before any BMRB 26544/26545 shift difference; descriptive', '']
    save = {}
    i51 = SITE - 1
    nm = np.array([k - 1 for k in NMOTIF])
    for T in ('278K', '300K'):
        xw, iw = load(T, 'WT')
        xp, ip = load(T, 'HIPK1T51')
        Cw, Cp = contacts(xw), contacts(xp)
        pw, pp = Cw.mean(0), Cp.mean(0)
        dp = pp - pw
        r = np.abs(dp).sum(1)
        net = dp.sum(1)
        h = len(xw) // 2
        noise = 0.5 * (np.abs(Cw[:h].mean(0) - Cw[h:].mean(0)).sum(1) + np.abs(Cp[:h].mean(0) - Cp[h:].mean(0)).sum(1))
        save.update({f'{T}_pw': pw, f'{T}_pp': pp, f'{T}_r': r, f'{T}_net': net, f'{T}_noise_halves': noise})
        L.append(f'== {T}: {iw["n_frames_total"]} frames per run, {len(xw)} production frames used')
        L.append(f'  Rg: WT {rg(xw).mean():.3f} nm, pT51 {rg(xp).mean():.3f} nm ({(rg(xp).mean() / rg(xw).mean() - 1) * 100:+.1f}%)')
        L.append(f'  total contacts per frame: WT {pw.sum() / 2:.1f}, pT51 {pp.sum() / 2:.1f}')
        L.append(f'  contacts of residue 51 with the N-motif (4-12), summed frequency: WT {pw[i51, nm].sum():.3f}, '
                 f'pT51 {pp[i51, nm].sum():.3f}')
        acid = np.array([k for k, c in enumerate(SEQ) if c in 'DE' and k not in nm])
        L.append(f'  N-motif contacts with all D/E residues, summed frequency: WT {pw[np.ix_(nm, acid)].sum():.3f}, '
                 f'pT51 {pp[np.ix_(nm, acid)].sum():.3f}')
        L.append(f'  rewiring r_i = sum_j |dp_ij|: median {np.median(r):.3f}; within-run half-vs-half floor median '
                 f'{np.median(noise):.3f} (lower bound)')
        far = np.abs(np.arange(len(r)) - i51) > 3
        top = [k for k in np.argsort(-r) if far[k]][:10]
        L.append('  most rewired residues more than 3 from T51: ' + ', '.join(f'{SEQ[k]}{k + 1} ({r[k]:.2f}, net {net[k]:+.2f})' for k in top))
        iu = np.triu_indices(len(r), 3)
        order = np.argsort(-np.abs(dp[iu]))[:8]
        L.append('  largest single contact changes: ' + ', '.join(
            f'{SEQ[iu[0][k]]}{iu[0][k] + 1}-{SEQ[iu[1][k]]}{iu[1][k] + 1} {dp[iu][k]:+.3f}' for k in order))
        reg = {'N-motif 4-12': nm, '13-40': np.arange(12, 40), '41-47': np.arange(40, 47), '55-75': np.arange(54, 75),
               '76-102': np.arange(75, 102)}
        L.append('  mean r_i by region: ' + ', '.join(f'{k} {r[v].mean():.3f}' for k, v in reg.items()))
        L.append('')
        if T == '300K':
            phi, _ = pn.page4_phi(xw, SEQ)
            x = -pn.sp.BETA_C * pn.DQ * phi[:, i51]
            w = np.exp(x - x.max())
            w /= w.sum()
            prw = np.tensordot(w, Cw.astype(float), axes=1)
            rrw = np.abs(prw - pw).sum(1)
            netrw = (prw - pw).sum(1)
            save.update({'300K_r_rw': rrw, '300K_net_rw': netrw, '300K_w': w})
            L.append(f'  reweighting (WT run, dq -2 at T51): effective sample fraction {1 / (w ** 2).sum() / len(w):.3f}; '
                     f'Rg {(w @ rg(xw) / rg(xw).mean() - 1) * 100:+.1f}%')
            L.append(f'  reweighted vs direct, per residue: Spearman r_i {spearmanr(rrw, r).correlation:+.3f}, '
                     f'net contact change {spearmanr(netrw, net).correlation:+.3f}')
            L.append(f'  reweighted contacts of residue 51 with the N-motif: {prw[i51, nm].sum():.3f} (WT {pw[i51, nm].sum():.3f})')
            L.append('')
    np.savez_compressed(os.path.join(OUT, 'page4_prediction_existing_runs.npz'), seq=SEQ, **save)
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'page4_prediction_existing_runs.txt'), 'w').write(txt)
    print(txt)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    k = np.arange(1, len(SEQ) + 1)
    for T, c in (('278K', '#1f77b4'), ('300K', '#d62728')):
        ax[0].plot(k, save[f'{T}_r'], color=c, lw=1.4, label=f'direct, {T}')
        ax[1].plot(k, save[f'{T}_net'], color=c, lw=1.4, label=f'direct, {T}')
    ax[0].plot(k, save['300K_noise_halves'], color='0.6', lw=1, ls=':', label='half-vs-half floor, 300K')
    ax[0].plot(k, save['300K_r_rw'], color='k', lw=1, ls='--', label='reweighted, 300K')
    ax[1].plot(k, save['300K_net_rw'], color='k', lw=1, ls='--', label='reweighted, 300K')
    for a in ax:
        a.axvspan(3.5, 12.5, color='#2ca02c', alpha=0.12, lw=0)
        a.axvspan(47.5, 54.5, color='0.85', lw=0)
        a.axvline(SITE, color='0.3', lw=0.8)
    ax[1].axhline(0, color='0.5', lw=0.6)
    ax[0].set_ylabel('contact rewiring  sum |dp|')
    ax[1].set_ylabel('net contact change  sum dp')
    ax[1].set_xlabel('PAGE4 residue  (green: basic N-motif 4-12; grey: T51 +-3, excluded in the NMR test)')
    ax[0].legend(fontsize=8, frameon=False, ncol=2)
    ax[0].set_title('PAGE4 T51 phosphorylation, CALVADOS-2 prediction from BENDER-MT runs (before the NMR comparison)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'page4_prediction_existing_runs.png'), dpi=150)


if __name__ == '__main__':
    main()
