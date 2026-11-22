"""
Do salt, pH and phosphorylation retune the ensemble? Responses from existing 300 K trajectories with the same identity
as for temperature: d<A>/d(lambda) = -beta Cov(A, dH/d(lambda)), A = Rg.

Salt (ionic strength I): H depends on I only through the screening constant kappa, kappa ~ sqrt(I), so
  dH/d ln I = (kappa / 2) dH/dkappa.  Reported as % change of <Rg> per doubling of I (x ln 2).
  CALVADOS-2: dV_DH/dkappa = sum q_i q_j eps (e^(-kappa rc) - e^(-kappa r)) for r < rc = 4 nm. In part A it is
  recovered per frame from the stored FDT energies, dV_DH/dkappa = (dVdT_DH - (eps'/eps) V_DH) / kappa'.
  Mpipi-GG: dE_DH/dkappa = -sum 332.06 q_i q_j e^(-kappa r) / 80 for r < the Coulomb cutoff (kcal/mol A).
pH (CALVADOS-2 only; Mpipi-GG charges do not depend on pH): only His titrates, q = 1/(1 + 10^(pH - 6)) = 0.031 at
  pH 7.5, so dH/dpH = (dq/dpH) sum_His phi_h with dq/dpH = -ln10 q (1 - q); phi_i = dV_DH/dq_i is the electrostatic
  potential at residue i from all other charges. Reported as % change of <Rg> per pH unit.
Phosphorylation (charge model, every Ser and Thr): adding dq at residue i changes the energy by exactly dq * phi_i
  per frame (V_DH is bilinear in the charges), so both the linear response -beta Cov(Rg, dq phi_i) and the exact
  reweighting <Rg w>/<w>, w = exp(-beta dq phi_i), are available; the Kish effective sample size of w says when
  the reweighting can be trusted. dq = -2 in CALVADOS-2 and -1.5 in Mpipi-GG (its charges are scaled by 0.75).
  Charge only: a phosphate's size and stickiness are not modelled.

Part A (--part npz): CALVADOS-2, every chain, 300 K, salt from the npz energies (no trajectories).
Part B (--part traj): the proteins simulated in both models, 300 K: salt (both models; CALVADOS-2 also checks part A),
  pH (CALVADOS-2), phosphorylation per Ser/Thr (both models).
Output: Traj Level Analysis/Analysis_FDT/salt_ph_phospho/
Usage: python3 salt_ph_phospho.py --part npz|traj|summary [--workers 8] [--limit N]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import TODO_T2_trajectory_loader as L  # noqa: E402
import TODO_T2_fdt_computation as fdt  # noqa: E402
import mpipi_temperature_hpc as mh  # noqa: E402
import basin_susceptibility as bs  # noqa: E402

bc = L.bc
T = 300.0
OUT_DIR = os.path.join(L.OUT_ROOT, 'salt_ph_phospho')
CAL = fdt.get_dielectric_params(T)
BETA_C = 1.0 / (fdt.R_GAS * T)
BETA_M = 1.0 / (bs.R_KCAL * T)
LN2 = np.log(2.0)
DQ = {'CALVADOS-2': -2.0, 'Mpipi-GG': -1.5}
N_BLOCKS = 10
MAX_ELEMS = 4_000_000


def cov(a, b):
    return float(((a - a.mean()) * (b - b.mean())).mean())


def jk_se(a, b, f):
    """Delete-one-block jackknife SE of f(a, b)."""
    blocks = np.array_split(np.arange(len(a)), N_BLOCKS)
    v = np.array([f(np.delete(a, bl), np.delete(b, bl)) for bl in blocks])
    return float(np.sqrt((N_BLOCKS - 1) / N_BLOCKS * ((v - v.mean()) ** 2).sum()))


def salt_slope(rg, dHdlnI, beta):
    """% change of <Rg> per doubling of I, and its block-jackknife SE."""
    f = lambda a, b: -beta * cov(a, b) / a.mean() * 100 * LN2  # noqa: E731
    return f(rg, dHdlnI), jk_se(rg, dHdlnI, f)


# ---------------------------------------------------------------- part A: CALVADOS-2 salt from the npz energies

def npz_salt(uid):
    p = os.path.join(L.OUT_ROOT, 'npz', '300K', f'{uid}.npz')
    try:
        z = np.load(p)
        V, dVdT, rg = (z[k].astype(float) for k in ('V_DH', 'dVdT_DH', 'rg'))
        dVdk = (dVdT - CAL['deps_yu_dT'] / CAL['eps_yu'] * V) / CAL['dkappa_dT']
        s, se = salt_slope(rg, CAL['kappa'] / 2 * dVdk, BETA_C)
        return {'UniProt_ID': uid, 'salt_rg_pct_per_doubling': s, 'salt_se': se, 'V_DH_mean': float(V.mean())}
    except Exception as e:  # noqa: BLE001
        return {'UniProt_ID': uid, 'error': repr(e)}


def part_npz(workers, limit):
    ids = pd.read_csv(os.path.join(L.OUT_ROOT, 'fdt_results_300K.csv'), usecols=['UniProt_ID', 'metric']).query("metric == 'rg'").UniProt_ID.tolist()
    ids = ids[:limit] if limit else ids
    print(f'part A: {len(ids)} CALVADOS-2 chains from npz', flush=True)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers * 2) as ex:
        rows = list(ex.map(npz_salt, ids))
    A = pd.DataFrame(rows)
    A.to_csv(os.path.join(OUT_DIR, 'salt_calvados_all.csv'), index=False)
    ok = A.dropna(subset=['salt_rg_pct_per_doubling'])
    print(f'  done: {len(ok)} ok, {A.error.notna().sum() if "error" in A else 0} failed; median {ok.salt_rg_pct_per_doubling.median():.3f} % per doubling', flush=True)


# ---------------------------------------------------------------- part B: trajectories, both models

def site_matrix(i, j, q, n):
    """Sparse (n x pairs): phi = W @ M.T gives each residue's potential from per-pair weights W (no self term)."""
    p = np.arange(len(i))
    return sparse.csr_matrix((np.concatenate([q[j], q[i]]), (np.concatenate([i, j]), np.concatenate([p, p]))), shape=(n, len(i)))


def calvados_terms(uid, hint):
    f = bs.frames_calvados(uid, 300, hint)
    if f is None:
        raise FileNotFoundError('no 300 K CALVADOS-2 frames')
    xyz, _, rg = f
    z = np.load(os.path.join(L.OUT_ROOT, 'npz', '300K', f'{uid}.npz'))
    seq = str(z['sequence'])
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
    i, j, *_ = fdt.pair_setup(q, lam, sig)
    M = site_matrix(i, j, q, len(seq))
    k, eps, rc = CAL['kappa'], CAL['eps_yu'], fdt.RC_YU
    ec = np.exp(-k * rc)
    nf = len(xyz)
    V, dVdk, phi = np.zeros(nf), np.zeros(nf), np.zeros((nf, len(seq)))
    step = max(1, MAX_ELEMS // max(1, len(i)))
    for a in range(0, nf, step):
        r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
        inside = r < rc
        er = np.exp(-k * r)
        Y = np.where(inside, er / r - ec / rc, 0.0)
        V[a:a + step] = (q[i] * q[j] * eps * Y).sum(1)
        dVdk[a:a + step] = (q[i] * q[j] * eps * np.where(inside, ec - er, 0.0)).sum(1)
        phi[a:a + step] = (M @ (eps * Y).T).T
    vz = z['V_DH'].astype(float)
    if np.corrcoef(V, vz)[0, 1] < 0.999 or abs(V.mean() - vz.mean()) > 1e-3 * max(1.0, abs(vz.mean())) + 1e-3:
        raise ValueError('recomputed V_DH does not match the npz')
    dVdk_npz = (z['dVdT_DH'].astype(float) - CAL['deps_yu_dT'] / CAL['eps_yu'] * vz) / CAL['dkappa_dT']
    return seq, rg, CAL['kappa'] / 2 * dVdk, CAL['kappa'] / 2 * dVdk_npz, phi, BETA_C


def mpipi_terms(uid):
    f = bs.frames_mpipi(uid, 300)
    if f is None:
        raise FileNotFoundError('no 300 K Mpipi-GG frames')
    xyz, _, rg = f
    z = np.load(os.path.join(bs.MP_DIR, 'npz_300K', f'{uid}.npz'))
    seq = str(z['sequence'])
    sysd = mh.read_mpipi_system(os.path.join(bs.MP_RUNS[300], uid))
    pp = mh.mpipi_pairs(sysd)
    i, j = pp['i'], pp['j']
    q = np.array([sysd['qtype'].get(int(t), 0.0) for t in sysd['types']])
    rc_all = max(sysd['coul'].values())
    rc_c = pp['rc_c']
    M = site_matrix(i, j, q, len(q))
    k, pref = sysd['kappa'], mh.QQR2E / sysd['diel']
    nf = len(xyz)
    E, dEdk, phi = np.zeros(nf), np.zeros(nf), np.zeros((nf, len(q)))
    step = max(1, MAX_ELEMS // max(1, len(i)))
    for a in range(0, nf, step):
        r = 10.0 * np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)       # A
        er = np.exp(-k * r)
        E[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * er / r, 0.0).sum(1)
        dEdk[a:a + step] = np.where(r < rc_c, -pref * pp['qq'] * er, 0.0).sum(1)
        phi[a:a + step] = (M @ np.where(r < rc_all, pref * er / r, 0.0).T).T
    ez = z['E_DH'].astype(float)
    if np.corrcoef(E, ez)[0, 1] < 0.999 and np.abs(ez).max() > 1e-6:
        raise ValueError('recomputed E_DH does not match the npz')
    return seq, rg, k / 2 * dEdk, None, phi, BETA_M


def analyse(task):
    uid, model, hint = task
    try:
        seq, rg, dHdlnI, dHdlnI_npz, phi, beta = calvados_terms(uid, hint) if model == 'CALVADOS-2' else mpipi_terms(uid)
        s, se = salt_slope(rg, dHdlnI, beta)
        row = {'UniProt_ID': uid, 'model': model, 'N': len(seq), 'rg_mean': float(rg.mean()),
               'salt_rg_pct_per_doubling': s, 'salt_se': se,
               'FCR': sum(seq.count(a) for a in 'DEKR') / len(seq),
               'NCPR': (seq.count('K') + seq.count('R') - seq.count('D') - seq.count('E')) / len(seq),
               'n_His': seq.count('H'), 'n_ST': seq.count('S') + seq.count('T')}
        if dHdlnI_npz is not None:
            row['salt_from_npz'] = salt_slope(rg, dHdlnI_npz, beta)[0]
        if model == 'CALVADOS-2':
            his = [p for p, c in enumerate(seq) if c == 'H']
            if his:
                qh = fdt.his_charge(fdt.PH)
                dHdpH = -np.log(10) * qh * (1 - qh) * phi[:, his].sum(1)
                f = lambda a, b: -beta * cov(a, b) / a.mean() * 100  # noqa: E731
                row['pH_rg_pct_per_unit'] = f(rg, dHdpH)
                row['pH_se'] = jk_se(rg, dHdpH, f)
        sites = []
        dq = DQ[model]
        for p, c in enumerate(seq):
            if c not in 'ST':
                continue
            dH = dq * phi[:, p]
            lin = -beta * cov(rg, dH) / rg.mean() * 100
            x = -beta * dH
            w = np.exp(x - x.max())
            rw = (float((w * rg).sum() / w.sum()) / rg.mean() - 1) * 100
            ess = float(w.sum() ** 2 / (w ** 2).sum() / len(w))
            sites.append({'UniProt_ID': uid, 'model': model, 'pos': p + 1, 'residue': c, 'phi_mean': float(phi[:, p].mean()),
                          'dRg_pct_linear': lin, 'dRg_pct_reweight': rw, 'ess_frac': ess})
        if sites:
            S = pd.DataFrame(sites)
            row['phospho_max_abs_pct'] = float(S.dRg_pct_reweight.abs().max())
            row['phospho_median_abs_pct'] = float(S.dRg_pct_reweight.abs().median())
        return row, sites, None
    except Exception as e:  # noqa: BLE001
        return None, None, f'{model} {uid}: {e!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def part_traj(workers, limit):
    per = pd.read_csv(bs.PER, index_col=0).dropna(subset=['chi_C', 'chi_M', 'nu_C', 'nu_M'])
    kingdom = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = [(u, m, kingdom.get(u)) for m in ('CALVADOS-2', 'Mpipi-GG') for u in per.index]
    tasks = tasks[:limit] if limit else tasks
    print(f'part B: {len(tasks)} chains, {workers} workers', flush=True)
    rows, sites, fails = [], [], []
    with mp.Pool(workers) as pool:
        for n, (r, s, err) in enumerate(pool.imap_unordered(analyse, tasks), 1):
            if err:
                fails.append(err)
            else:
                rows.append(r)
                sites += s
            if n % 50 == 0:
                print(f'  {n}/{len(tasks)} done, {len(fails)} failed', flush=True)
    P = pd.DataFrame(rows).join(per[['nu_C', 'nu_M']], on='UniProt_ID')
    P.to_csv(os.path.join(OUT_DIR, 'per_protein_common.csv'), index=False)
    pd.DataFrame(sites).to_csv(os.path.join(OUT_DIR, 'phospho_sites_common.csv'), index=False)
    with open(os.path.join(OUT_DIR, 'failed_common.txt'), 'w') as fh:
        fh.write('\n'.join(fails))
    print(f'  done: {len(P)} ok, {len(fails)} failed', flush=True)
    for f in fails[:5]:
        print('   ', f[:300])


# ---------------------------------------------------------------- summary of parts A and B

def charge_class(fcr, ncpr):
    """Das-Pappu regions: weak polyampholytes, strong polyampholytes, polyelectrolytes (and the boundary zone)."""
    if abs(ncpr) > 0.35:
        return 'polyelectrolyte'
    if fcr < 0.25:
        return 'weak polyampholyte'
    if fcr > 0.35:
        return 'strong polyampholyte'
    return 'boundary'


def part_summary(workers=None, limit=None):
    from scipy import stats
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 30)
    out = []
    A = pd.read_csv(os.path.join(OUT_DIR, 'salt_calvados_all.csv'))
    bio = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'sequence', 'SCD', 'kappa'], low_memory=False).set_index('UniProt_ID')
    per = pd.read_csv(bs.PER, index_col=0)
    A = A.join(bio, on='UniProt_ID').join(per[['chi_C', 'nu_C']], on='UniProt_ID')
    seq = A.sequence
    A['FCR'] = seq.map(lambda q: sum(q.count(a) for a in 'DEKR') / len(q))
    A['NCPR'] = seq.map(lambda q: (q.count('K') + q.count('R') - q.count('D') - q.count('E')) / len(q))
    A['cls'] = [charge_class(f, n) for f, n in zip(A.FCR, A.NCPR)]
    A['sig'] = (A.salt_rg_pct_per_doubling.abs() > 2 * A.salt_se)
    s = A.salt_rg_pct_per_doubling
    out.append(('A all CALVADOS-2: n', len(A)))
    out.append(('A median salt % per doubling', s.median()))
    out.append(('A share expanding / compacting (significant, |s| > 2 SE)', f"{(A.sig & (s > 0)).mean():.3f} / {(A.sig & (s < 0)).mean():.3f}"))
    out.append(('A 5th / 95th percentile', f'{s.quantile(0.05):.2f} / {s.quantile(0.95):.2f}'))
    for c in ('NCPR', 'FCR', 'SCD', 'kappa', 'nu_C', 'chi_C'):
        out.append((f'A Spearman salt vs {c}', stats.spearmanr(s, A[c], nan_policy='omit')[0]))
    out.append(('A Spearman salt vs |NCPR|', stats.spearmanr(s, A.NCPR.abs())[0]))
    g = A.groupby('cls').salt_rg_pct_per_doubling
    cls = pd.DataFrame({'n': g.size(), 'median': g.median(), 'q25': g.quantile(0.25), 'q75': g.quantile(0.75),
                        'share_expand': A.groupby('cls').apply(lambda x: (x.sig & (x.salt_rg_pct_per_doubling > 0)).mean(), include_groups=False),
                        'share_compact': A.groupby('cls').apply(lambda x: (x.sig & (x.salt_rg_pct_per_doubling < 0)).mean(), include_groups=False)})
    cls.to_csv(os.path.join(OUT_DIR, 'salt_by_charge_class.csv'))
    B = pd.read_csv(os.path.join(OUT_DIR, 'per_protein_common.csv'))
    bc_ = B[B.model == 'CALVADOS-2'].set_index('UniProt_ID')
    bm = B[B.model == 'Mpipi-GG'].set_index('UniProt_ID')
    j = bc_[['salt_rg_pct_per_doubling', 'salt_se', 'salt_from_npz', 'NCPR', 'FCR']].join(bm[['salt_rg_pct_per_doubling', 'salt_se']], rsuffix='_M')
    out.append(('B npz vs direct salt (CALVADOS-2), max abs diff', float((j.salt_rg_pct_per_doubling - j.salt_from_npz).abs().max())))
    out.append(('B salt Spearman CALVADOS-2 vs Mpipi-GG (common set)', stats.spearmanr(j.salt_rg_pct_per_doubling, j.salt_rg_pct_per_doubling_M)[0]))
    both = (j.salt_rg_pct_per_doubling.abs() > 2 * j.salt_se) & (j.salt_rg_pct_per_doubling_M.abs() > 2 * j.salt_se_M)
    out.append(('B same sign where both significant', f'{(np.sign(j.salt_rg_pct_per_doubling) == np.sign(j.salt_rg_pct_per_doubling_M))[both].mean():.3f} of {both.sum()}'))
    out.append(('B median salt CALVADOS-2 / Mpipi-GG', f'{j.salt_rg_pct_per_doubling.median():.3f} / {j.salt_rg_pct_per_doubling_M.median():.3f}'))
    ph = bc_.dropna(subset=['pH_rg_pct_per_unit'])
    out.append(('B pH (CALVADOS-2, chains with His): n', len(ph)))
    out.append(('B pH median |% per unit| / max |%|', f'{ph.pH_rg_pct_per_unit.abs().median():.3f} / {ph.pH_rg_pct_per_unit.abs().max():.3f}'))
    out.append(('B pH share significant (|x| > 2 SE)', float((ph.pH_rg_pct_per_unit.abs() > 2 * ph.pH_se).mean())))
    out.append(('B pH Spearman |pH effect| vs n_His/N', stats.spearmanr(ph.pH_rg_pct_per_unit.abs(), ph.n_His / ph.N)[0]))
    S = pd.read_csv(os.path.join(OUT_DIR, 'phospho_sites_common.csv'))
    for m, x in S.groupby('model'):
        good = x[x.ess_frac >= 0.1]
        out.append((f'B phospho {m}: sites / with ESS >= 0.1', f'{len(x)} / {len(good)}'))
        out.append((f'B phospho {m}: median |dRg| % (reweight)', good.dRg_pct_reweight.abs().median()))
        out.append((f'B phospho {m}: share |dRg| > 2% / > 5%', f'{(good.dRg_pct_reweight.abs() > 2).mean():.3f} / {(good.dRg_pct_reweight.abs() > 5).mean():.3f}'))
        out.append((f'B phospho {m}: share expanding (dRg > 0)', float((good.dRg_pct_reweight > 0).mean())))
        out.append((f'B phospho {m}: linear vs reweight Pearson', float(np.corrcoef(good.dRg_pct_linear, good.dRg_pct_reweight)[0, 1])))
        out.append((f'B phospho {m}: Spearman dRg vs chain NCPR', stats.spearmanr(good.dRg_pct_reweight, good.join(B[B.model == m].set_index('UniProt_ID').NCPR, on='UniProt_ID').NCPR)[0]))
        out.append((f'B phospho {m}: Spearman dRg vs local potential phi', stats.spearmanr(good.dRg_pct_reweight, good.phi_mean)[0]))
    sc = S[S.model == 'CALVADOS-2'].set_index(['UniProt_ID', 'pos'])
    sm = S[S.model == 'Mpipi-GG'].set_index(['UniProt_ID', 'pos'])
    ss = sc[['dRg_pct_reweight', 'ess_frac']].join(sm[['dRg_pct_reweight', 'ess_frac']], rsuffix='_M', how='inner')
    ss = ss[(ss.ess_frac >= 0.1) & (ss.ess_frac_M >= 0.1)]
    out.append(('B phospho site agreement across models: Spearman (sites)', f'{stats.spearmanr(ss.dRg_pct_reweight, ss.dRg_pct_reweight_M)[0]:.3f} (n {len(ss)})'))
    out.append(('B phospho same sign across models', float((np.sign(ss.dRg_pct_reweight) == np.sign(ss.dRg_pct_reweight_M)).mean())))
    top = ss.assign(mean=(ss.dRg_pct_reweight + ss.dRg_pct_reweight_M) / 2).reindex(columns=['dRg_pct_reweight', 'dRg_pct_reweight_M', 'mean'])
    top = top.reindex(top['mean'].abs().sort_values(ascending=False).index).head(15)
    top.to_csv(os.path.join(OUT_DIR, 'phospho_top_sites_common.csv'))
    T = pd.DataFrame(out, columns=['quantity', 'value'])
    T.to_csv(os.path.join(OUT_DIR, 'summary.csv'), index=False)
    print(T.to_string(index=False))
    print(cls.round(3).to_string())
    print(top.round(2).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['npz', 'traj', 'summary'], required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    {'npz': part_npz, 'traj': part_traj, 'summary': part_summary}[a.part](a.workers, a.limit)
    if a.part == 'summary':
        sys.exit(0)
    bc.write_provenance(os.path.join(OUT_DIR, 'salt_calvados_all.csv' if a.part == 'npz' else 'per_protein_common.csv'), __file__,
                        {'T': T, 'ionic_M': CAL['ionic'], 'kappa_cal_nm': CAL['kappa'], 'pH': fdt.PH, 'his_pKa': fdt.HIS_PKA,
                         'phospho_dq': DQ, 'part': a.part, 'blocks': N_BLOCKS}, [bs.PER])
