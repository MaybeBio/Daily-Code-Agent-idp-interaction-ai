"""
TODO-T2: trajectory loader and runner for the BENDER-MT FDT susceptibility
==========================================================================
Connects the BENDER-MT CALVADOS-2 DCDs to TODO_T2_fdt_computation.py.

FRAMES (identical to the reanalysis, bender_common.load_calvados)
  production = second half of each DCD (frames[n//2:]), whatever its length
  K_FRAMES = 1000 evenly strided production frames -> energies, G_T and every observable
  (--graph-frames N evaluates the graph metrics on N of them instead; `check` uses 50, the frames
   behind the network means in Analysis_<T>K_v2, to confirm the definitions reproduce them)

OBSERVABLES (definitions of Traj Level Analysis/Analysis_<T>K_v2)
  rg, ree, delta   per frame as in bender_common.geometry (unweighted; delta = 1.5 sum(ev-mean)^2 / tr^2)
  graph metrics    Calpha-Calpha < 0.8 nm with |i-j| >= 3 (bender_common.adjacency) and the
                   graph_features.py definitions: fragmentation_index = |LCC|/N, global_efficiency,
                   avg_clustering, transitivity, degree_assortativity (NaN when undefined)

PARAMETERS  calvados_residue_params() from residues_CALVADOS2.csv, using the topology's sequence.

OUTPUT (Traj Level Analysis/Analysis_FDT/)
  fdt_results_<T>K.csv   one row per (protein, observable): chi_total (+ block SE), both
                         decompositions, <A> (+ block SE) from the same frames
  npz/<T>K/<ID>.npz      per-frame energies, G_T and observables (float32), for bootstraps and pair maps
                         (graph_pos: positions of the graph frames among frame_index)
  fdt_validation.csv     all five temperatures: chi_T and local slope d<A>/dT per T, integration check
  fdt_validation_summary.csv
  logs/                  failed_<T>K.tsv (and the run log when started as below)

USAGE
  python3 TODO_T2_trajectory_loader.py check [--temp 310] [--uid A0A014QYM5]
  python3 TODO_T2_trajectory_loader.py run --temp 278 [--workers 4] [--ids IDS.txt] [--limit N] [--redo]
  python3 TODO_T2_trajectory_loader.py all [--workers 4]      # 278, 288, 310, 320, 300 K, then validate
  python3 TODO_T2_trajectory_loader.py validate
Reruns (e.g. the 411 truncated 300 K trajectories once ingested):
  python3 TODO_T2_trajectory_loader.py run --temp 300 --ids reruns_300K_ids.txt --redo
"""

import sys

sys.dont_write_bytecode = True  # keep the NAS script folder free of __pycache__

import argparse  # noqa: E402
import hashlib  # noqa: E402
import os  # noqa: E402
import platform  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
import warnings  # noqa: E402
import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import networkx as nx  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
REANALYSIS_DIR = os.environ.get('BENDER_REANALYSIS_DIR', '/Users/taseefrahman/Downloads/Scripts/bender_reanalysis')
sys.path.insert(0, HERE)
sys.path.insert(0, REANALYSIS_DIR)
import TODO_T2_fdt_computation as fdt  # noqa: E402
import bender_common as bc  # noqa: E402

if platform.system() == 'Darwin':
    mp.set_start_method('fork', force=True)

# ===========================================================================
# CONFIGURATION
# ===========================================================================

BENDER_MT = '/Volumes/User Homes/Nu Project/BENDER-MT'
TRAJ_ROOT_FMT = BENDER_MT + '/BENDER_{T}K'
V2_DIR_FMT = BENDER_MT + '/Traj Level Analysis/Analysis_{T}K_v2'
OUT_ROOT = BENDER_MT + '/Traj Level Analysis/Analysis_FDT'
BIO_CSV = bc.BIO_CSV
RESIDUES_CSV_PATH = '/Users/taseefrahman/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'
if not os.path.exists(RESIDUES_CSV_PATH):   # this Mac's home folder (2026-10-10)
    RESIDUES_CSV_PATH = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')

TEMPS = [278, 288, 300, 310, 320]
RUN_ORDER = [278, 288, 310, 320, 300]  # 300 K last: organise_release.py may still refile BENDER_300K
K_FRAMES = 1000
GRAPH_FRAMES = 1000   # all K frames; contacts decorrelate within one saved frame, so these are ~independent
GEOMETRY = ['rg', 'ree', 'delta']
GRAPH_METRICS = ['fragmentation_index', 'global_efficiency', 'avg_clustering', 'transitivity',
                 'degree_assortativity']
ENERGY_KEYS = ['V_AH', 'V_DH', 'V_bond', 'dVdT_DH', 'dudT_DH', 'G_T']


def load_residue_params(csv_path: str = RESIDUES_CSV_PATH) -> pd.DataFrame:
    """residues_CALVADOS2.csv indexed by one-letter code (columns lambdas, sigmas, q, bondlength)."""
    return pd.read_csv(csv_path).set_index('one')


# ===========================================================================
# TRAJECTORIES AND PER-FRAME OBSERVABLES
# ===========================================================================

def is_flat_layout(base: str, ids: list) -> bool:
    """BENDER_<T>/<ID>/ (278, 288, 310, 320 K) rather than BENDER_300K/<group>/<ID>/ (as in the v2 runner)."""
    return any(os.path.isdir(os.path.join(base, u)) for u in ids[:50])


def find_traj_dir(uid: str, base: str, hint: str = None, flat: bool = None):
    """
    Folder holding <uid>.dcd and top.pdb: BENDER_<T>/<uid> (flat) or BENDER_300K/<group>/<uid>.
    Group folders are listed when the protein is processed, so a refiled 300 K folder is still found:
    the BENDER_BIO kingdom first, 'Other' last (it duplicates the lineage folders). In a flat layout
    a missing ID is reported at once (never scan ~10k protein folders as if they were groups).
    """
    has = lambda d: os.path.isfile(os.path.join(d, f'{uid}.dcd')) and os.path.isfile(os.path.join(d, 'top.pdb'))  # noqa: E731
    if has(os.path.join(base, uid)):
        return os.path.join(base, uid)
    if flat is None:
        flat = is_flat_layout(base, [uid])
    if flat:
        return None
    groups = sorted(g for g in os.listdir(base) if os.path.isdir(os.path.join(base, g)))
    groups = ([hint] if hint in groups else []) + [g for g in groups if g not in (hint, 'Other')] + \
             (['Other'] if 'Other' in groups else [])
    for g in groups:
        if has(os.path.join(base, g, uid)):
            return os.path.join(base, g, uid)
    return None


def geometry_per_frame(xyz: np.ndarray) -> dict:
    """Per-frame Rg, Ree and asphericity with the bender_common.geometry definitions."""
    r = xyz - xyz.mean(axis=1, keepdims=True)
    ev = np.linalg.eigvalsh(np.einsum('fni,fnj->fij', r, r) / xyz.shape[1])
    tr = ev.sum(axis=1)
    return {
        'rg': np.sqrt(tr),
        'ree': np.linalg.norm(xyz[:, -1] - xyz[:, 0], axis=1),
        'delta': 1.5 * ((ev - ev.mean(axis=1, keepdims=True)) ** 2).sum(axis=1) / (tr ** 2 + 1e-12),
    }


def graph_per_frame(xyz: np.ndarray) -> dict:
    """Per-frame contact-graph metrics with the graph_features.py definitions (|i-j| >= 3 graphs)."""
    n_frames, n_res, _ = xyz.shape
    out = {m: np.full(n_frames, np.nan) for m in GRAPH_METRICS}
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


def analyse(uid: str, d: str, T: int, residues_df: pd.DataFrame, npz_dir: str,
            k_frames: int = K_FRAMES, graph_frames: int = GRAPH_FRAMES, timings: dict = None) -> list:
    """All FDT rows for one protein at one temperature; writes npz/<T>K/<uid>.npz."""
    t = [time.time()]
    xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=k_frames)
    t.append(time.time())
    q, lam, sig, b0 = fdt.calvados_residue_params(seq, residues_df)
    E = fdt.per_frame_energies(xyz, q, lam, sig, b0, fdt.get_dielectric_params(float(T)))
    t.append(time.time())
    geo = geometry_per_frame(xyz)
    gi = bc.even_indices(len(xyz), graph_frames)
    graph = graph_per_frame(xyz[gi])
    t.append(time.time())
    if timings is not None:
        timings.update(load=t[1] - t[0], energies=t[2] - t[1], observables=t[3] - t[2])

    base = {'UniProt_ID': uid, 'T': T, 'n_residues': len(seq), 'n_frames_dcd': info['n_frames_total'],
            'production_start': info['production_start'], 'stride': info['stride'],
            'max_bond_nm': bc.max_bond(xyz)}
    rows = [{**base, 'metric': m, **fdt.susceptibility(a, E, T)} for m, a in geo.items()]
    Eg = {k: v[gi] for k, v in E.items()}
    rows += [{**base, 'metric': m, **fdt.susceptibility(a, Eg, T)} for m, a in graph.items()]

    f32 = lambda x: np.asarray(x, dtype=np.float32)  # noqa: E731
    np.savez_compressed(
        os.path.join(npz_dir, f'{uid}.npz'), T=T, sequence=seq,
        frame_index=info['production_start'] + info['stride'] * np.arange(len(xyz)), graph_pos=gi,
        **{k: f32(E[k]) for k in ENERGY_KEYS}, **{m: f32(a) for m, a in geo.items()},
        **{f'graph_{m}': f32(a) for m, a in graph.items()})
    return rows


def _worker(args):
    uid, base, flat, T, npz_dir, hint, k_frames, graph_frames = args
    try:
        d = find_traj_dir(uid, base, hint, flat)
        if d is None:
            return uid, None, 'missing trajectory'
        return uid, analyse(uid, d, T, RESIDUES, npz_dir, k_frames, graph_frames), None
    except Exception as e:
        return uid, None, f'{e!r} | {traceback.format_exc()[-400:]}'.replace('\n', ' ')


RESIDUES = None  # set in run() before the pool forks


# ===========================================================================
# RUNNER
# ===========================================================================

def _append_rows(csv_path: str, rows: list, tries: int = 6):
    """Single writer (the parent): append rows keeping the header's column order. Transient SMB
    errors are retried with back-off, as in traj_analysis_calvados_v2._append."""
    df = pd.DataFrame(rows)
    for k in range(tries):
        try:
            if os.path.exists(csv_path):
                cols = pd.read_csv(csv_path, nrows=0).columns.tolist()
                if [c for c in df.columns if c not in cols]:
                    pd.concat([pd.read_csv(csv_path), df], ignore_index=True).to_csv(csv_path, index=False)
                else:
                    df.reindex(columns=cols).to_csv(csv_path, mode='a', header=False, index=False)
            else:
                df.to_csv(csv_path, index=False)
            return
        except OSError as e:
            if k == tries - 1:
                raise
            print(f'  append retry {k + 1} after {e!r}', flush=True)
            time.sleep(5 * (k + 1))


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def run(T: int, out: str = OUT_ROOT, ids_file: str = None, workers: int = 4, limit: int = None,
        redo: bool = False, k_frames: int = K_FRAMES, graph_frames: int = GRAPH_FRAMES,
        base_dir: str = None):
    global RESIDUES
    RESIDUES = load_residue_params()
    base = base_dir or TRAJ_ROOT_FMT.format(T=T)
    npz_dir = os.path.join(out, 'npz', f'{T}K')
    os.makedirs(npz_dir, exist_ok=True)
    os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    csv_path = os.path.join(out, f'fdt_results_{T}K.csv')
    fail_path = os.path.join(out, 'logs', f'failed_{T}K.tsv')

    bio = pd.read_csv(BIO_CSV, usecols=['UniProt_ID', 'kingdom'])
    bio['UniProt_ID'] = bio['UniProt_ID'].astype(str)
    hint = dict(zip(bio['UniProt_ID'], bio['kingdom']))
    ids = ([l.strip() for l in open(ids_file) if l.strip()] if ids_file else bio['UniProt_ID'].tolist())

    done = set()
    if os.path.exists(csv_path):
        prev = pd.read_csv(csv_path)
        if redo:
            keep = prev[~prev['UniProt_ID'].astype(str).isin(ids)]
            keep.to_csv(csv_path, index=False)
            print(f'--redo: removed {len(prev) - len(keep)} rows for {len(ids)} requested IDs', flush=True)
            prev = keep
        done = set(prev['UniProt_ID'].astype(str))
    todo = [u for u in ids if u not in done][:limit]
    flat = is_flat_layout(base, ids)
    print(f'{T} K: {len(ids)} requested, {len(done)} done, {len(todo)} to run with {workers} workers '
          f'(K={k_frames}, graph frames={graph_frames}) from {base} ({"flat" if flat else "group folders"})',
          flush=True)

    jobs = [(u, base, flat, T, npz_dir, hint.get(u), k_frames, graph_frames) for u in todo]
    n_fail = 0
    t0 = time.time()
    with mp.Pool(workers) as pool:
        for i, (uid, rows, err) in enumerate(pool.imap_unordered(_worker, jobs), 1):
            if rows:
                _append_rows(csv_path, rows)
            else:
                n_fail += 1
                with open(fail_path, 'a') as fh:
                    fh.write(f'{uid}\t{err}\n')
                print(f'  FAILED {uid}: {err[:160]}', flush=True)
            if i % 100 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'[{T} K {i}/{len(jobs)}] {el / 60:.1f} min, ~{el / i * (len(jobs) - i) / 3600:.1f} h left, '
                      f'{n_fail} failed', flush=True)

    bc.write_provenance(csv_path, __file__, {
        'T': T, 'K_FRAMES': k_frames, 'GRAPH_FRAMES': graph_frames, 'production': 'frames[n//2:]',
        'contact_cutoff_nm': bc.CONTACT_CUTOFF_NM, 'contact_min_sep': bc.CONTACT_MIN_SEP,
        'fdt_module': os.path.join(HERE, 'TODO_T2_fdt_computation.py'),
        'fdt_module_sha256': _sha256(os.path.join(HERE, 'TODO_T2_fdt_computation.py')),
        'residues_csv_sha256': _sha256(RESIDUES_CSV_PATH), 'workers': workers, 'base_dir': base,
        'eps_ah': fdt.EPS_AH, 'rc_ah': fdt.RC_AH, 'rc_yu': fdt.RC_YU, 'kb_bond': fdt.KB_BOND,
        'ionic': fdt.IONIC_STRENGTH, 'pH': fdt.PH, 'n_blocks': fdt.N_BLOCKS_CHI,
    }, [BIO_CSV, RESIDUES_CSV_PATH, base] + ([ids_file] if ids_file else []))
    print(f'{T} K finished: {n_fail} failed (see {fail_path})' if n_fail else f'{T} K finished', flush=True)


# ===========================================================================
# VALIDATION ACROSS THE FIVE TEMPERATURES
# ===========================================================================

def load_results(out: str = OUT_ROOT, columns: list = None) -> pd.DataFrame:
    """The five fdt_results_<T>K.csv tables stacked, with `truncated` = DCD shorter than a complete run
    (bc.expected_frames; all temperatures used the same setup, some 278 K runs are longer)."""
    tabs = []
    for T in TEMPS:
        path = os.path.join(out, f'fdt_results_{T}K.csv')
        if not os.path.exists(path):
            raise FileNotFoundError(f'{path} missing; run --temp {T} first')
        cols = None if columns is None else list(dict.fromkeys(
            ['UniProt_ID', 'T', 'metric', 'n_residues', 'n_frames_dcd'] + columns))
        tabs.append(pd.read_csv(path, usecols=cols))
    df = pd.concat(tabs, ignore_index=True).drop_duplicates(['UniProt_ID', 'metric', 'T'], keep='last')
    df['UniProt_ID'] = df['UniProt_ID'].astype(str)
    df['truncated'] = df['n_frames_dcd'] < df['n_residues'].map(bc.expected_frames)
    return df


def validate(out: str = OUT_ROOT):
    """
    Per protein and observable: chi_T and the local slope d<A>/dT at each temperature, and the
    integration check. Every sequence with all five temperatures is kept; truncated trajectories
    are flagged, not excluded.
    """
    temps = np.array(TEMPS, dtype=float)
    df = load_results(out, ['chi_total', 'chi_total_se', 'A_mean', 'A_mean_se'])
    wide = df.pivot_table(index=['UniProt_ID', 'metric'], columns='T',
                          values=['chi_total', 'chi_total_se', 'A_mean', 'A_mean_se'])
    wide = wide.dropna(subset=[(v, T) for v in ('chi_total', 'A_mean') for T in TEMPS])

    trunc = {T: set(df.loc[(df['T'] == T) & df['truncated'], 'UniProt_ID']) for T in TEMPS}
    rows = []
    for (uid, metric), r in wide.iterrows():
        chi = np.array([r[('chi_total', T)] for T in TEMPS])
        chi_se = np.array([r[('chi_total_se', T)] for T in TEMPS])
        A = np.array([r[('A_mean', T)] for T in TEMPS])
        A_se = np.array([r[('A_mean_se', T)] for T in TEMPS])
        row = {'UniProt_ID': uid, 'metric': metric}
        for T, c, s in zip(TEMPS, chi, fdt.local_slopes(A, temps)):
            row[f'chi_{T}K'] = c
            row[f'slope_{T}K'] = s
        row.update(fdt.global_integration_check(chi, chi_se, A, A_se, temps))
        row['linear_slope'] = fdt.linear_slope(A, temps)
        row['truncated_temps'] = ';'.join(str(T) for T in TEMPS if uid in trunc[T])
        rows.append(row)
    val = pd.DataFrame(rows)
    path = os.path.join(out, 'fdt_validation.csv')
    val.to_csv(path, index=False)

    summ = []
    for metric, g in val.groupby('metric'):
        for T in TEMPS:
            x, y = g[f'chi_{T}K'], g[f'slope_{T}K']
            d = y - x
            summ.append({'metric': metric, 'T': T, 'n': len(g), 'mean_chi': x.mean(), 'mean_slope': y.mean(),
                         'mean_slope_minus_chi': d.mean(), 'se_slope_minus_chi': d.std() / np.sqrt(len(d)),
                         'pearson_r': x.corr(y), 'spearman_r': x.corr(y, method='spearman')})
        summ.append({'metric': metric, 'T': 'integral', 'n': len(g), 'mean_chi': g['integral'].mean(),
                     'mean_slope': g['observed_diff'].mean(),
                     'mean_slope_minus_chi': (g['observed_diff'] - g['integral']).mean(),
                     'se_slope_minus_chi': (g['observed_diff'] - g['integral']).std() / np.sqrt(len(g)),
                     'pearson_r': g['integral'].corr(g['observed_diff']),
                     'spearman_r': g['integral'].corr(g['observed_diff'], method='spearman'),
                     'frac_abs_z_lt_2': (g['z'].abs() < 2).mean()})
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(out, 'fdt_validation_summary.csv'), index=False)
    bc.write_provenance(path, __file__, {'temps': TEMPS, 'slopes': 'np.gradient on the T grid',
                                         'integral': 'trapezoid, SEs in quadrature'},
                        [os.path.join(out, f'fdt_results_{T}K.csv') for T in TEMPS])
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(summ.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print(f'saved {path} ({len(val)} rows)')


# ===========================================================================
# CHECK: sanity checks + agreement with the reanalysis on one protein
# ===========================================================================

def check(T: int = 310, uid: str = 'A0A014QYM5'):
    fdt.run_sanity_checks()
    residues = load_residue_params()
    base = TRAJ_ROOT_FMT.format(T=T)
    d = find_traj_dir(uid, base)
    timings = {}
    with tempfile.TemporaryDirectory() as tmp:  # 50 graph frames = the reanalysis's network frames
        rows = pd.DataFrame(analyse(uid, d, T, residues, tmp, graph_frames=50, timings=timings)).set_index('metric')
    print(f'{uid} at {T} K: timings (s) ' + ', '.join(f'{k} {v:.1f}' for k, v in timings.items()))
    v2 = pd.read_csv(os.path.join(V2_DIR_FMT.format(T=T), 'scalar_features_v2.csv')).set_index('UniProt_ID')
    ok = True
    for m in GEOMETRY + GRAPH_METRICS:
        if m in v2.columns and pd.notna(v2.loc[uid, m]):
            diff = abs(rows.loc[m, 'A_mean'] - v2.loc[uid, m])
            ok &= diff < 1e-9
            print(f'  <{m}> here {rows.loc[m, "A_mean"]:.10g} vs Analysis_{T}K_v2 {v2.loc[uid, m]:.10g}'
                  f'  {"✓" if diff < 1e-9 else "✗"}')
    cols = ['chi_total', 'chi_total_se', 'chi_Boltzmann', 'chi_dielectric', 'f_dielectric', 'chi_AH',
            'chi_bond', 'chi_DH', 'n_frames']
    with pd.option_context('display.width', 200):
        print(rows[cols].to_string(float_format=lambda v: f'{v:.3e}'))
    print('check passed' if ok else 'check FAILED: observables differ from the reanalysis')
    return ok


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['check', 'run', 'all', 'validate'])
    ap.add_argument('--temp', type=int, choices=TEMPS)
    ap.add_argument('--uid', default='A0A014QYM5')
    ap.add_argument('--out', default=OUT_ROOT)
    ap.add_argument('--ids', default=None, help='file with one UniProt ID per line (default: BENDER_BIO.csv)')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--redo', action='store_true', help='recompute IDs already in the results table')
    ap.add_argument('--k-frames', type=int, default=K_FRAMES)
    ap.add_argument('--graph-frames', type=int, default=GRAPH_FRAMES)
    ap.add_argument('--base-dir', default=None, help='trajectory root (default BENDER-MT/BENDER_<T>K)')
    a = ap.parse_args()
    if a.mode == 'check':
        sys.exit(0 if check(a.temp or 310, a.uid) else 1)
    elif a.mode == 'run':
        if a.temp is None:
            ap.error('run needs --temp')
        run(a.temp, a.out, a.ids, a.workers, a.limit, a.redo, a.k_frames, a.graph_frames, a.base_dir)
    elif a.mode == 'all':
        for T in RUN_ORDER:
            run(T, a.out, a.ids, a.workers, a.limit, a.redo, a.k_frames, a.graph_frames)
        validate(a.out)
    else:
        validate(a.out)
