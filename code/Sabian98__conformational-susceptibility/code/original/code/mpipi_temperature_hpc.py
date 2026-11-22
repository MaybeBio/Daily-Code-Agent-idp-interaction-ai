"""
mpipi_temperature_hpc.py — analyse the MPIPI-GG runs of the 992 paired sequences at one temperature,
on GAIVI, with the same conventions as the 300 K MPIPI-GG analysis (mpipi_analysis_v2.py +
mpipi_graph_v2.py), so that only a small results table has to be copied back.

Per run directory (seq.dat, traj.xtc, run.in):
  * the thermostat temperature is read from run.in and checked against --temp
  * all traj.xtc frames are production; coordinates are unwrapped from minimum-image bond vectors
  * geometry (Rg, Ree, asphericity), primary nu (RMS per pair -> mean per s -> fit s >= 5) with R^2,
    curve-fit SE and 5-block SE, and nu on each half of the run: all frames
  * contact map (C-alpha < 0.8 nm, |i-j| >= 3) on 1,000 evenly spaced frames; pi-pi and cation-pi
    contact frequencies (mean per-pair frequency over F/W/Y pairs and Arg-F/W/Y pairs, as in
    Scripts/dynamic_pi_score.py) with pair counts
  * contact-network features (graph_features.py) on 50 of those frames
  * FDT susceptibility on the same 1,000 frames: per-frame Rg, Ree, asphericity and the five network metrics
    (definitions of the CALVADOS-2 FDT runner) and the per-frame potential energy from log.lammps (thermo `pe`,
    written at the xtc steps). Mpipi-GG's Hamiltonian does not depend on temperature (fixed dielectric 80 and
    Debye kappa), so chi_T(A) = d<A>/dT = Cov(A, PotEng) / (R T^2) exactly; block-jackknife SE over 10 blocks.
    Energy-term split chi_WF + chi_DH + chi_bond (Derivation Eq. 11.14) from the Wang-Frenkel, Debye-Hueckel and
    bond energies recomputed per frame from run.in/seq.dat (Joseph et al. 2021, Eqs. 1-5); their sum is checked
    against PotEng per run (terms_minus_pe_rms_kcal, terms_pe_r). chi on each half of the frames (equilibration).

Outputs: --out (ensemble means, one row per run), <out>_fdt.csv (one row per run and observable, the columns of
BENDER-MT's fdt_results_<T>K.csv) and <out dir>/npz_<T>K/<ID>.npz (per-frame energies and observables).

Copy three files to one folder on GAIVI: this script, bender_common.py and graph_features.py (the last two from
Scripts/bender_reanalysis; numpy, scipy, pandas, mdtraj, networkx required). Then submit each temperature with
submit_mpipi_temperature.sh (same folder), or run directly:
  python /path/to/mpipi_temperature_hpc.py --runs /path/to/runs_278K --temp 278 --out /data/taseefr/.../mpipi_278K.csv --workers 16
and copy the CSVs (and .provenance.json) back to BENDER-BIO/MPIPI/MPIPI Results/, the npz folder to
BENDER-MT/Traj Level Analysis/Analysis_FDT_MPIPI/.
"""

import argparse
import multiprocessing as mp
import os
import re
import sys
import time
import traceback

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True  # the script folder may be read-only (HPC) or the NAS code folder
# shared modules: next to this script on GAIVI (copied flat), else bender_reanalysis on the Mac
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', '/Users/taseefrahman/Downloads/Scripts/bender_reanalysis'))
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import bender_common as bc  # noqa: E402
from graph_features import graph_features_from_adj  # noqa: E402

K_FRAMES, GRAPH_FRAMES = 1000, 50          # as traj_analysis_calvados_v2 / mpipi_graph_v2
AROMATIC, CATION_PI = set('FWY'), set('R')  # as dynamic_pi_score.py
NETWORK = ['fragmentation_index', 'global_efficiency', 'avg_clustering', 'transitivity', 'degree_assortativity']
R_KCAL = 1.98720425864083e-3               # kcal mol^-1 K^-1 (LAMMPS units real)
N_BLOCKS_CHI = 10                          # as the CALVADOS-2 FDT


def run_temperature(run_in):
    m = re.search(r'fix\s+\S+\s+all\s+langevin\s+([\d.]+)\s+([\d.]+)', open(run_in).read())
    return float(m.group(1)) if m else float('nan')


def xtc_steps(xtc):
    import mdtraj as md
    with md.formats.XTCTrajectoryFile(xtc) as fh:
        return fh.read(atom_indices=[0])[2].astype(np.int64)


def thermo_potential_energy(log):
    """step -> PotEng (kcal/mol) from every thermo block of log.lammps that reports PotEng (later blocks win)."""
    pe, cols = {}, None
    for line in open(log):
        s = line.split()
        if s and s[0] == 'Step':
            cols = s if 'PotEng' in s else None
            continue
        if cols and len(s) == len(cols) and s[0].isdigit():
            try:
                pe[int(s[0])] = float(s[cols.index('PotEng')])
            except ValueError:
                pass
        elif cols and line.startswith('Loop time'):
            cols = None
    return pe


def geometry_per_frame(xyz):
    """Per-frame Rg, Ree, asphericity with the bender_common.geometry definitions."""
    r = xyz - xyz.mean(axis=1, keepdims=True)
    ev = np.linalg.eigvalsh(np.einsum('fni,fnj->fij', r, r) / xyz.shape[1])
    tr = ev.sum(axis=1)
    return {'rg': np.sqrt(tr), 'ree': np.linalg.norm(xyz[:, -1] - xyz[:, 0], axis=1),
            'delta': 1.5 * ((ev - ev.mean(axis=1, keepdims=True)) ** 2).sum(axis=1) / (tr ** 2 + 1e-12)}


def graph_per_frame(xyz):
    """Per-frame network metrics with the graph_features.py definitions (|i-j| >= 3 graphs)."""
    import warnings
    import networkx as nx
    n_frames, n_res, _ = xyz.shape
    out = {m: np.full(n_frames, np.nan) for m in NETWORK}
    for f in range(n_frames):
        G = nx.from_numpy_array(bc.adjacency(xyz[f:f + 1])[0].astype(int))
        out['fragmentation_index'][f] = len(max(nx.connected_components(G), key=len)) / n_res
        out['global_efficiency'][f] = nx.global_efficiency(G)
        out['avg_clustering'][f] = nx.average_clustering(G)
        out['transitivity'][f] = nx.transitivity(G)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                out['degree_assortativity'][f] = nx.degree_assortativity_coefficient(G)
        except Exception:
            pass
    return out


# ---------------------------------------------------------------------------
# Mpipi-GG energy terms per frame (Joseph et al. 2021 Methods Eqs. 1-5; parameters from each run's run.in/seq.dat)
#   bonds  K (r - r0)^2 (LAMMPS harmonic, K = k/2);  Debye-Hueckel  C q_i q_j exp(-kappa r) / (eps_r r), r < rc;
#   Wang-Frenkel  eps alpha [(sigma/r)^2mu - 1] [(rc/r)^2mu - 1]^2nu, r < rc;  special_bonds fene: 1-2 pairs excluded
# The three terms reproduce LAMMPS's thermo PotEng frame by frame up to xtc rounding (checked, r > 0.999).
# ---------------------------------------------------------------------------
QQR2E = 332.06371       # kcal A / (mol e^2), LAMMPS real units
MAX_ELEMS = 4e6         # frames x pairs per chunk


def read_mpipi_system(d):
    txt = open(f'{d}/run.in').read()
    style = re.search(r'^pair_style.*$', txt, re.M).group(0)
    sysd = {'diel': float(re.search(r'^dielectric\s+(\S+)', txt, re.M).group(1)),
            # sequences without charged residues are run with pair_style wf/cut only (no Debye-Hueckel term)
            'kappa': float(re.search(r'coul/debye\s+(\S+)', style).group(1)) if 'coul/debye' in style else 0.0,
            'qtype': {int(t): float(q) for t, q in re.findall(r'^set type\s+(\d+)\s+charge\s+(\S+)', txt, re.M)},
            # 'wf/cut' keyword only under pair_style hybrid/overlay (charged sequences)
            'wf': {(int(i), int(j)): (float(e), float(s), int(nu), int(mu), float(rc)) for i, j, e, s, nu, mu, rc in
                   re.findall(r'^pair_coeff\s+(\d+)\s+(\d+)\s+(?:wf/cut\s+)?([-+.\deE]+)\s+([-+.\deE]+)\s+(\d+)\s+(\d+)'
                              r'\s+([-+.\deE]+)\s*$', txt, re.M)},
            'coul': {(int(i), int(j)): float(rc) for i, j, rc in
                     re.findall(r'^pair_coeff\s+(\d+)\s+(\d+)\s+coul/debye\s+(\S+)', txt, re.M)}}
    sysd['kb'], sysd['r0'] = map(float, re.search(r'^bond_coeff\s+1\s+(\S+)\s+(\S+)', txt, re.M).groups())
    if 'special_bonds fene' not in txt:
        raise ValueError('expected special_bonds fene (1-2 pairs excluded)')
    lines = open(f'{d}/seq.dat').read().splitlines()
    ia = next(k for k, l in enumerate(lines) if l.strip() == 'Atoms')
    ib = next(k for k, l in enumerate(lines) if l.strip() == 'Bonds')
    atoms = sorted((int(s[0]), int(s[2])) for s in (l.split() for l in lines[ia + 1:ib]) if len(s) >= 7)
    sysd['types'] = np.array([t for _, t in atoms])
    sysd['bonds'] = np.array([(int(s[2]) - 1, int(s[3]) - 1) for s in (l.split() for l in lines[ib + 1:]) if len(s) >= 4])
    return sysd


def mpipi_pairs(sysd):
    t = sysd['types']
    i, j = np.triu_indices(len(t), k=1)
    bonded = {(min(a, b), max(a, b)) for a, b in sysd['bonds']}
    keep = np.array([(a, b) not in bonded for a, b in zip(i, j)])
    i, j = i[keep], j[keep]
    key = lambda a, b: (min(a, b), max(a, b))  # noqa: E731
    P = np.array([sysd['wf'][key(t[a], t[b])] for a, b in zip(i, j)], dtype=float)
    q = np.array([sysd['qtype'].get(tt, 0.0) for tt in t])
    rc_c = np.array([sysd['coul'].get(key(t[a], t[b]), 0.0) for a, b in zip(i, j)])
    eps, sig, nu, mu, rc = P.T
    alpha = 2 * nu * (rc / sig) ** (2 * mu) * ((1 + 2 * nu) / (2 * nu * ((rc / sig) ** (2 * mu) - 1))) ** (2 * nu + 1)
    return {'i': i, 'j': j, 'eps_alpha': eps * alpha, 'sig': sig, 'nu': nu, 'mu': mu, 'rc': rc,
            'qq': q[i] * q[j], 'rc_c': rc_c}


def mpipi_energy_terms(xyz_nm, box_nm, sysd, pp):
    """Per-frame E_WF, E_DH, E_bond (kcal/mol) with minimum-image distances."""
    x = np.asarray(xyz_nm, dtype=float) * 10.0
    L = float(box_nm) * 10.0
    n = len(x)
    E = {k: np.zeros(n) for k in ('E_WF', 'E_DH', 'E_bond')}
    step = max(1, int(MAX_ELEMS // max(1, len(pp['i']))))
    for a in range(0, n, step):
        dv = x[a:a + step, pp['i']] - x[a:a + step, pp['j']]
        dv -= L * np.round(dv / L)
        r = np.linalg.norm(dv, axis=2)
        wf = pp['eps_alpha'] * ((pp['sig'] / r) ** (2 * pp['mu']) - 1) * ((pp['rc'] / r) ** (2 * pp['mu']) - 1) ** (2 * pp['nu'])
        E['E_WF'][a:a + step] = np.where(r < pp['rc'], wf, 0.0).sum(1)
        dh = QQR2E * pp['qq'] * np.exp(-sysd['kappa'] * r) / (sysd['diel'] * r)
        E['E_DH'][a:a + step] = np.where(r < pp['rc_c'], dh, 0.0).sum(1)
    bv = x[:, sysd['bonds'][:, 0]] - x[:, sysd['bonds'][:, 1]]
    bv -= L * np.round(bv / L)
    E['E_bond'] = (sysd['kb'] * (np.linalg.norm(bv, axis=2) - sysd['r0']) ** 2).sum(1)
    return E


def _cov(x, y):
    return float(((x - x.mean()) * (y - y.mean())).sum() / (len(x) - 1))


def susceptibility(a, pe, temp, terms=None, n_blocks=N_BLOCKS_CHI):
    """
    chi = Cov(A, PotEng)/(R T^2) (exact: H does not depend on T) with a delete-one-block jackknife SE, the
    energy-term split chi_WF + chi_DH + chi_bond (Derivation Eq. 11.14; equals chi_total up to xtc rounding of the
    recomputed terms), chi on each half of the frames (equilibration check), and <A> with its block SE.
    Frames with undefined A are dropped. chi_Boltzmann = chi_total and chi_dielectric = 0 by construction.
    """
    ok = np.isfinite(a)
    a, pe = a[ok], pe[ok]
    terms = {k: v[ok] for k, v in (terms or {}).items()}
    n = len(a)
    keys = ['chi_total', 'chi_total_se', 'chi_Boltzmann', 'chi_dielectric', 'chi_WF', 'chi_DH', 'chi_bond',
            'chi_half1', 'chi_half2', 'A_mean', 'A_mean_se', 'A_sd']
    if n < 2 * n_blocks:
        return {**{k: np.nan for k in keys}, 'n_frames': n}
    rt2 = R_KCAL * temp ** 2
    chi = lambda i: _cov(a[i], pe[i]) / rt2  # noqa: E731
    blocks = np.array_split(np.arange(n), n_blocks)
    jk = np.array([chi(np.setdiff1d(np.arange(n), b)) for b in blocks])
    bm = np.array([a[b].mean() for b in blocks])
    out = {'chi_total': chi(np.arange(n)),
           'chi_total_se': float(np.sqrt((n_blocks - 1) / n_blocks * ((jk - jk.mean()) ** 2).sum()))}
    out['chi_Boltzmann'], out['chi_dielectric'] = out['chi_total'], 0.0
    for k, name in (('E_WF', 'chi_WF'), ('E_DH', 'chi_DH'), ('E_bond', 'chi_bond')):
        out[name] = _cov(a, terms[k]) / rt2 if k in terms else np.nan
    h = n // 2
    out['chi_half1'], out['chi_half2'] = chi(np.arange(h)), chi(np.arange(h, n))
    out.update({'A_mean': float(a.mean()), 'A_mean_se': float(bm.std(ddof=1) / np.sqrt(n_blocks)),
                'A_sd': float(a.std(ddof=1)), 'n_frames': n})
    return out


def fdt_rows(uid, d, seq, xyz, box_nm, temp, npz_dir):
    """FDT rows on the K_FRAMES contact frames; writes <npz_dir>/<uid>.npz."""
    steps = xtc_steps(f'{d}/traj.xtc')
    pe_map = thermo_potential_energy(f'{d}/log.lammps')
    kidx = bc.even_indices(len(xyz), K_FRAMES)
    pe = np.array([pe_map.get(int(s), np.nan) for s in steps[kidx]])
    matched = np.isfinite(pe)
    if matched.mean() < 0.99:
        raise ValueError(f'only {matched.sum()}/{len(pe)} frames have a thermo PotEng at their step')
    xk = xyz[kidx][matched]
    pe = pe[matched]
    sysd = read_mpipi_system(d)
    if len(sysd['types']) != xyz.shape[1]:
        raise ValueError(f'seq.dat has {len(sysd["types"])} atoms, xtc {xyz.shape[1]}')
    E = mpipi_energy_terms(xk, box_nm, sysd, mpipi_pairs(sysd))
    resid = E['E_WF'] + E['E_DH'] + E['E_bond'] - pe
    obs = {**geometry_per_frame(xk), **graph_per_frame(xk)}
    base = {'UniProt_ID': uid, 'T': temp, 'n_residues': len(seq), 'n_frames_xtc': len(xyz),
            'pe_mean_kcal': float(pe.mean()), 'pe_sd_kcal': float(pe.std(ddof=1)),
            'terms_minus_pe_rms_kcal': float(np.sqrt((resid ** 2).mean())),
            'terms_pe_r': float(np.corrcoef(pe, pe + resid)[0, 1])}
    rows = [{**base, 'metric': m, **susceptibility(a, pe, temp, E)} for m, a in obs.items()]
    if npz_dir:
        f32 = lambda x: np.asarray(x, dtype=np.float32)  # noqa: E731
        np.savez_compressed(os.path.join(npz_dir, f'{uid}.npz'), T=temp, sequence=seq,
                            step=steps[kidx][matched], PotEng_kcal=f32(pe), **{k: f32(v) for k, v in E.items()},
                            **{m: f32(a) for m, a in obs.items()})
    return rows


def pair_mean(F, idx_a, idx_b):
    pairs = {(min(i, j), max(i, j)) for i in idx_a for j in idx_b if i != j and abs(i - j) >= bc.CONTACT_MIN_SEP}
    if not pairs:
        return 0.0, 0
    ii, jj = np.array(sorted(pairs)).T
    return float(F[ii, jj].mean()), len(pairs)


def analyse(args):
    d, temp, npz_dir = args      # passed per task: workers do not inherit globals under spawn/forkserver
    uid = os.path.basename(d.rstrip('/'))
    t0 = time.time()
    fdt = []
    try:
        seq = bc.read_mpipi_sequence(f'{d}/seq.dat')
        t_run = run_temperature(f'{d}/run.in')
        xyz, info = bc.load_mpipi(f'{d}/traj.xtc', n_res=len(seq))          # all frames, unwrapped
        if info['n_frames_total'] < 100:
            return {'UniProt_ID': uid, 'status': f"too_few_frames:{info['n_frames_total']}"}
        row = {'UniProt_ID': uid, 'temperature_K': temp, 'run_in_temperature_K': t_run,
               'n_residues': len(seq), 'n_frames': info['n_frames_total'], 'box_nm': info['box_nm'],
               'max_bond_nm': round(info['max_bond_nm'], 4)}
        row.update(bc.geometry(xyz))
        f = bc.flory(xyz)
        row.update({k: f[k] for k in ('nu', 'A0', 'nu_fit_err', 'nu_r2', 'nu_block_se', 'nu_legacy', 'nu_legacy_r2')})
        h = len(xyz) // 2
        row['nu_half1'] = bc.nu_primary_from_xyz(xyz[:h])['nu']
        row['nu_half2'] = bc.nu_primary_from_xyz(xyz[h:])['nu']
        xk = xyz[bc.even_indices(len(xyz), K_FRAMES)]
        F = bc.contact_frequency_map(xk)
        aro = [i for i, a in enumerate(seq) if a in AROMATIC]
        cat = [i for i, a in enumerate(seq) if a in CATION_PI]
        row['pi_pi_contact_freq'], row['n_aromatic_pairs'] = pair_mean(F, aro, aro) if len(aro) >= 2 else (0.0, 0)
        row['cation_pi_contact_freq'], row['n_cation_pi_pairs'] = pair_mean(F, cat, aro) if cat and aro else (0.0, 0)
        xg = xk[bc.even_indices(len(xk), GRAPH_FRAMES)]
        g, _ = graph_features_from_adj(bc.adjacency(xg), seq)
        row.update({k: g[k] for k in NETWORK})
        row.update({f'{k}_std': g[f'{k}_std'] for k in NETWORK})
        row['n_frames_graph'] = len(xg)
        row['sequence'] = seq
        try:
            fdt = fdt_rows(uid, d, seq, xyz, info['box_nm'], temp, npz_dir)
            row['fdt_status'] = 'ok'
        except Exception as e:
            row['fdt_status'] = f'failed: {e!r}'[:200]
        row['seconds'] = round(time.time() - t0, 1)
        row['status'] = 'ok' if abs(t_run - temp) < 0.5 else f'temperature_mismatch:{t_run}'
        return row, fdt
    except Exception as e:
        return {'UniProt_ID': uid, 'status': f'failed: {e!r}'[:300], 'trace': traceback.format_exc()[-400:]}, fdt


def main(runs, temp, out, workers, limit=None, ids=None):
    temp = float(temp)
    names = [l.strip() for l in open(ids) if l.strip()] if ids else os.listdir(runs)
    dirs = sorted(os.path.join(runs, u) for u in names
                  if os.path.exists(os.path.join(runs, u, 'traj.xtc')) and os.path.exists(os.path.join(runs, u, 'seq.dat')))
    dirs = dirs[:limit]
    out = os.path.abspath(out)
    out_fdt = out[:-4] + '_fdt.csv' if out.endswith('.csv') else out + '_fdt.csv'
    npz_dir = os.path.join(os.path.dirname(out), f'npz_{temp:g}K')
    os.makedirs(npz_dir, exist_ok=True)
    print(f'{len(dirs)} runs at {temp:g} K, {workers} workers', flush=True)
    t0, rows, fdt = time.time(), [], []
    with mp.Pool(workers) as pool:
        for i, (r, f) in enumerate(pool.imap_unordered(analyse, [(d, temp, npz_dir) for d in dirs]), 1):
            rows.append(r)
            fdt += f
            if i % 50 == 0 or r.get('status') != 'ok' or r.get('fdt_status', 'ok') != 'ok':
                print(f"[{i}/{len(dirs)}] {r.get('status')}/fdt {r.get('fdt_status')}:{r['UniProt_ID']}  "
                      f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    pd.DataFrame(fdt).to_csv(out_fdt, index=False)
    params = {'temperature_K': temp, 'frames': 'all (unwrapped) for geometry/nu; '
              f'{K_FRAMES} for contacts and the FDT; {GRAPH_FRAMES} for the graph means',
              'contact_cutoff_nm': bc.CONTACT_CUTOFF_NM, 'contact_min_sep': bc.CONTACT_MIN_SEP,
              'fdt': 'chi = Cov(A, thermo PotEng)/(R T^2), R = 1.98720425864083e-3 kcal/mol/K, '
                     f'{N_BLOCKS_CHI}-block jackknife SE', 'npz_dir': npz_dir, 'workers': workers}
    for path in (out, out_fdt):
        bc.write_provenance(path, __file__, params, [runs])
    ok = df[df['status'] == 'ok']
    print(f'wrote {out}: {len(ok)} ok of {len(df)}; median nu {ok["nu"].median():.3f}; '
          f'max bond after unwrapping {ok["max_bond_nm"].max():.3f} nm')
    print(f'wrote {out_fdt}: {len(fdt)} rows; FDT ok for {int((df.get("fdt_status") == "ok").sum())} runs')
    print(df['status'].value_counts().head(10).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', required=True, help='folder with one sub-folder per sequence (seq.dat, traj.xtc, run.in)')
    ap.add_argument('--temp', required=True, type=float, help='temperature of this set, in K')
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=None, help='first N runs only (testing)')
    ap.add_argument('--ids', default=None, help='file with one run ID per line: analyse only these (e.g. late arrivals)')
    a = ap.parse_args()
    main(a.runs, a.temp, a.out, a.workers, a.limit, a.ids)
