"""
mpipi_recharged_fdt.py — FDT susceptibility for the Mpipi-Recharged (OpenMpipi / OpenMM) BENDER runs, on GAIVI.

Runs come from openmpipi_scripts/run_openmpipi_idp.py: <root>/<T>K/<ID>/<ID>.dcd + top.pdb, the CALVADOS-2 layout
(first half of the frames = equilibration, same frame counts, cubic box 50 nm or N/5 + 20 nm above 300 residues,
150 mM). The run log and state.xml are deleted after each run, so every frame's energy is recomputed here.

Mpipi-Recharged (Tejedor et al., ACS Cent. Sci. 2025): pair-specific Wang-Frenkel stickers, pair-specific asymmetric
Yukawa electrostatics A_ij exp(-kappa r)/r with A_ij at a fixed dielectric of 80, and kappa from the salt and the
temperature-dependent water dielectric, so H depends on T through the screening. Instead of assuming how OpenMpipi
encodes that, the system is rebuilt with OpenMpipi's own get_mpipi_system at T and at T +/- DT, and the energies of
the same frames are differenced:

    dH/dT = [H(x; T + DT) - H(x; T - DT)] / (2 DT)            (central difference, fixed coordinates)
    G_T   = d(beta H)/dT = -H / (R T^2) + (dH/dT) / (R T)
    chi_T(A) = -Cov(A, G_T)                                   (TODO_T2_fdt_computation.susceptibility)

Per frame, the energy of every OpenMM force is evaluated separately (one force group each) and labelled
wf / yukawa / bond / other from its class and energy expression; the labels and expressions are written to the
provenance. Decomposition passed to the shared susceptibility() (same keys as CALVADOS-2):
    V_AH <- Wang-Frenkel (sticker) energy     -> chi_WF
    V_DH <- Yukawa (+ any 'other') energy     -> chi_DH = Boltzmann part + screening part
    V_bond <- bonds                           -> chi_bond
    dVdT_DH <- dH/dT of the whole system      -> chi_dielectric (here: the screening term)
The script checks that dH/dT comes from the Yukawa force alone (dHdT_non_yukawa_frac).

Frames and observables exactly as the CALVADOS-2 FDT runner (TODO_T2_trajectory_loader.py): 1,000 evenly spaced
production frames (bender_common.load_calvados), per-frame Rg / Ree / asphericity and the five contact-network
metrics (C-alpha < 0.8 nm, |i-j| >= 3), 10-block jackknife SEs, chi on each half of the frames.

Modes:
  check --runs ROOT --temp T [--name ID]   one protein: force list and expressions, what changes in the serialized
                                           system between T and T + DT (how T enters), kappa change vs the paper's
                                           formula, group energies vs total, DT/2 vs DT agreement, timings
  logcheck --runs ROOT --temp T [--name ID] an unfinished run that still has <ID>.log and seg_*.dcd: recomputed
                                           total energy vs the potential energy the simulation itself logged, at the
                                           logged steps that fall on saved frames (proves the rebuilt system is the
                                           simulated one)
  run --runs ROOT --temp T --out DIR [--ids FILE] [--workers N] [--platform CPU|CUDA] [--limit N]
                                           every <ID> under ROOT/<T>K with <ID>.dcd and top.pdb; resumable
Outputs (run): DIR/recharged_fdt_<T>K.csv (one row per protein and observable), DIR/npz_<T>K/<ID>.npz (per-frame
energies, dH/dT and observables), DIR/logs/failed_<T>K.tsv. Copy the CSV to
BENDER-MT/Traj Level Analysis/Analysis_FDT_RECHARGED/ and the npz folder next to it.

On GAIVI: put this script, TODO_T2_fdt_computation.py, TODO_T2_trajectory_loader.py, bender_common.py and
graph_features.py in one folder; use the environment that runs run_openmpipi_idp.py (OpenMM + OpenMpipi + mdtraj +
networkx). Submit with submit_mpipi_recharged_fdt.sh. The script never changes directory.
"""
import sys

sys.dont_write_bytecode = True  # keep the code folder free of __pycache__

import argparse  # noqa: E402
import difflib  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import TODO_T2_trajectory_loader as L  # noqa: E402  (geometry_per_frame, graph_per_frame, frame settings)

fdt, bc = L.fdt, L.bc
METRICS = L.GEOMETRY + L.GRAPH_METRICS
K_FRAMES = L.K_FRAMES
DT = 1.0            # K, central difference for dH/dT
IONIC_M = 0.15      # run_openmpipi_idp.py default (150 mM)
EQUIL_FRAC = 0.5    # first 50% of the DCD frames discarded as equilibration (run_openmpipi_idp.py, CALVADOS prepare.py)
LABELS = ('wf', 'yukawa', 'bond', 'other')


def box_nm(n_res: int) -> float:
    """Cubic box edge used by run_openmpipi_idp.py."""
    return float(int(n_res / 5) + 20 if n_res > 300 else 50)


def build_system(name: str, seq: str, T: float, ionic_M: float, xyz0: np.ndarray):
    """The run's own system: OpenMpipi get_mpipi_system with the run_openmpipi_idp.py arguments."""
    import openmm as mm
    from openmm import unit
    from OpenMpipi import IDP, get_mpipi_system
    L_box = box_nm(len(seq))
    top = IDP(name, seq).topology
    box = [mm.Vec3(L_box, 0, 0), mm.Vec3(0, L_box, 0), mm.Vec3(0, 0, L_box)] * unit.nanometer
    top.setPeriodicBoxVectors(box)
    system = get_mpipi_system(xyz0, top, {name: []}, T, ionic_M * 1000, CM_remover=False, periodic=True)
    system.setDefaultPeriodicBoxVectors(*box)
    return system


def label_force(f) -> str:
    cls = f.__class__.__name__
    if cls in ('HarmonicBondForce', 'CustomBondForce'):
        return 'bond'
    if cls == 'CustomNonbondedForce':
        return 'yukawa' if 'exp(' in f.getEnergyFunction().replace(' ', '').lower() else 'wf'
    if cls == 'CMMotionRemover':
        return 'none'
    return 'other'


def describe_forces(system) -> list:
    out = []
    for k, f in enumerate(system.getForces()):
        expr = f.getEnergyFunction() if hasattr(f, 'getEnergyFunction') else ''
        glob = {}
        if hasattr(f, 'getNumGlobalParameters'):
            glob = {f.getGlobalParameterName(g): f.getGlobalParameterDefaultValue(g) for g in range(f.getNumGlobalParameters())}
        cutoff = f.getCutoffDistance().value_in_unit_system(__import__('openmm').unit.md_unit_system) \
            if hasattr(f, 'getCutoffDistance') else None
        out.append({'index': k, 'class': f.__class__.__name__, 'label': label_force(f), 'expression': expr[:300],
                    'globals': glob, 'cutoff_nm': cutoff})
    return out


class Evaluator:
    """Per-force energies (kJ/mol) of given frames under one system."""

    def __init__(self, system, platform='CPU'):
        import openmm as mm
        self.mm = mm
        self.forces = system.getForces()
        if len(self.forces) > 32:
            raise ValueError('more than 32 forces: cannot give each its own force group')
        for k, f in enumerate(self.forces):
            f.setForceGroup(k)
        props = {'Threads': '1'} if platform == 'CPU' else ({'Precision': 'double'} if platform == 'CUDA' else {})
        self.ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(platform), props)
        self.labels = [label_force(f) for f in self.forces]

    def energies(self, xyz: np.ndarray) -> np.ndarray:
        from openmm import unit
        out = np.zeros((len(xyz), len(self.forces)))
        kj = unit.kilojoule_per_mole
        for t, x in enumerate(xyz):
            self.ctx.setPositions(x)
            for k in range(len(self.forces)):
                out[t, k] = self.ctx.getState(getEnergy=True, groups=1 << k).getPotentialEnergy().value_in_unit(kj)
        return out

    def total(self, x: np.ndarray) -> float:
        from openmm import unit
        self.ctx.setPositions(x)
        return self.ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def run_schedule(n_res: int) -> dict:
    """Step arithmetic of run_openmpipi_idp.py (frames, save interval, segments, StateDataReporter interval)."""
    t_ns = 200 if n_res <= 150 else int(200 * (n_res / 150) ** 1.5)
    frames = 10000 if t_ns <= 300 else 15000 if t_ns <= 600 else 20000
    n_save = int((t_ns * 1000 / frames) / 0.02)
    seg_steps = frames // 10 * n_save
    return {'frames': frames, 'n_save': n_save, 'seg_frames': frames // 10, 'seg_steps': seg_steps,
            'log_every': max(n_save, seg_steps // 10)}


def logged_frames(run_dir: str, uid: str) -> tuple:
    """(xyz nm [M,N,3], logged potential energy kJ/mol [M], steps [M], seq) for every logged step that is a saved
    frame of a finished segment; the last log entry of a step wins (a re-run segment overwrites its DCD)."""
    import json
    import mdtraj as md
    top = os.path.join(run_dir, 'top.pdb')
    seq = ''.join(bc.AA3_TO_1.get(r.name[:3].upper(), 'X') for r in md.load_topology(top).residues)
    sch = run_schedule(len(seq))
    done = json.load(open(os.path.join(run_dir, 'progress.json')))['segments_done']
    pe = {}
    for line in open(os.path.join(run_dir, f'{uid}.log')):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        pe[int(f[0])] = float(f[1])
    xyz, E, steps = [], [], []
    for s in range(done):
        traj = md.load(os.path.join(run_dir, f'seg_{s:02d}.dcd'), top=top)
        for S in sorted(k for k in pe if s * sch['seg_steps'] < k <= (s + 1) * sch['seg_steps']):
            k, rem = divmod(S - s * sch['seg_steps'], sch['n_save'])
            if rem == 0 and 1 <= k <= traj.n_frames:
                xyz.append(traj.xyz[k - 1].astype(np.float64))
                E.append(pe[S])
                steps.append(S)
    return np.array(xyz), np.array(E), np.array(steps), seq


def logcheck(root: str, T: float, name: str = None, ionic_M: float = IONIC_M, platform: str = 'CPU'):
    tdir = os.path.join(root, f'{int(T)}K')
    cands = [u for u in sorted(os.listdir(tdir)) if os.path.isfile(os.path.join(tdir, u, f'{u}.log'))
             and os.path.isfile(os.path.join(tdir, u, 'progress.json'))]
    for uid in ([name] if name else cands):
        xyz, E_log, steps, seq = logged_frames(os.path.join(tdir, uid), uid)
        if not len(xyz):
            print(f'{uid}: no logged step on a saved frame')
            continue
        ev = Evaluator(build_system(uid, seq, T, ionic_M, xyz[0]), platform)
        E = np.array([ev.total(x) for x in xyz])
        d = E - E_log
        print(f'{uid}: N {len(seq)}, {len(E)} logged frames; recomputed - logged: mean {d.mean():+.4f}, '
              f'max |d| {np.abs(d).max():.4f} kJ/mol (energies {E_log.min():.1f} to {E_log.max():.1f}); '
              f'r = {np.corrcoef(E, E_log)[0, 1]:.6f}')


def load_run(run_dir: str, uid: str):
    xyz, seq, info = bc.load_calvados(os.path.join(run_dir, f'{uid}.dcd'), os.path.join(run_dir, 'top.pdb'), k=K_FRAMES)
    if info['production_start'] != int(info['n_frames_total'] * EQUIL_FRAC):   # load_calvados keeps frames n//2 onwards
        raise ValueError(f"production starts at frame {info['production_start']} of {info['n_frames_total']}, "
                         f'expected a {EQUIL_FRAC:.0%} equilibration discard')
    info['max_bond_raw_nm'] = bc.max_bond(xyz)
    if info['max_bond_raw_nm'] > 1.0:                 # a chain split by the periodic box: rebuild it whole
        xyz = bc.unwrap_chain(xyz, box_nm(len(seq)))
    info['max_bond_nm'] = bc.max_bond(xyz)
    if 'X' in seq:
        raise ValueError(f'unknown residue name in top.pdb: {seq}')
    return xyz, seq, info


def energies_at(name, seq, xyz, T, ionic_M, platform, dT=DT):
    """Per-force energies at T, and dH/dT per force by central difference at fixed coordinates."""
    ev = {t: Evaluator(build_system(name, seq, t, ionic_M, xyz[0]), platform) for t in (T, T - dT, T + dT)}
    labels = ev[T].labels
    if any(ev[t].labels != labels for t in ev):
        raise ValueError('force list differs between temperatures')
    E = {t: e.energies(xyz) for t, e in ev.items()}
    return labels, E[T], (E[T + dT] - E[T - dT]) / (2 * dT)


def per_frame(labels, E0, dEdT, T):
    lab = np.array(labels)
    s = lambda X, which: X[:, np.isin(lab, which)].sum(axis=1) if np.isin(lab, which).any() else np.zeros(len(X))  # noqa: E731
    RT = fdt.R_GAS * T
    V_WF, V_Y, V_bond, V_other = (s(E0, [x]) for x in LABELS)
    H = E0.sum(axis=1)
    dHdT = dEdT.sum(axis=1)
    V_DH = V_Y + V_other
    return {'V_AH': V_WF, 'V_DH': V_DH, 'V_bond': V_bond, 'H': H, 'dVdT_DH': dHdT,
            'dudT_DH': -V_DH / (RT * T) + dHdT / RT, 'G_T': -H / (RT * T) + dHdT / RT,
            'V_other': V_other, 'dHdT_yukawa': s(dEdT, ['yukawa']),
            'dHdT_non_yukawa': s(dEdT, ['wf', 'bond', 'other', 'none'])}


def analyse(args):
    uid, run_dir, T, ionic_M, platform, npz_dir = args
    t0 = time.time()
    try:
        xyz, seq, info = load_run(run_dir, uid)
        labels, E0, dEdT = energies_at(uid, seq, xyz, T, ionic_M, platform)
        e = per_frame(labels, E0, dEdT, T)
        geo = L.geometry_per_frame(xyz)
        graph = L.graph_per_frame(xyz)
        obs = {**geo, **graph}
        n = len(xyz)
        h1, h2 = np.arange(n // 2), np.arange(n // 2, n)
        nonY = np.sqrt(np.mean(e['dHdT_non_yukawa'] ** 2)) / max(np.sqrt(np.mean(e['dVdT_DH'] ** 2)), 1e-12)
        base = {'UniProt_ID': uid, 'T': T, 'n_residues': len(seq), 'n_frames_dcd': info['n_frames_total'],
                'expected_frames': bc.expected_frames(len(seq)),
                'truncated': info['n_frames_total'] < bc.expected_frames(len(seq)),
                'production_start': info['production_start'], 'stride': info['stride'],
                'max_bond_raw_nm': round(info['max_bond_raw_nm'], 4), 'max_bond_nm': round(info['max_bond_nm'], 4),
                'box_nm': box_nm(len(seq)), 'ionic_M': ionic_M, 'force_labels': ';'.join(labels),
                'H_mean_kJ': float(e['H'].mean()), 'H_sd_kJ': float(e['H'].std(ddof=1)),
                'dHdT_mean_kJ_per_K': float(e['dVdT_DH'].mean()), 'dHdT_non_yukawa_frac': float(nonY),
                'V_other_rms_kJ': float(np.sqrt(np.mean(e['V_other'] ** 2)))}
        keep = ['V_AH', 'V_DH', 'V_bond', 'H', 'dVdT_DH', 'dudT_DH', 'G_T']
        rows = []
        for m in METRICS:
            r = fdt.susceptibility(obs[m], {k: e[k] for k in keep}, T)
            r = {('chi_WF' if k == 'chi_AH' else 'chi_screening' if k == 'chi_dielectric'
                  else 'f_screening' if k == 'f_dielectric' else k): v for k, v in r.items()}
            r['chi_half1'] = fdt.susceptibility(obs[m][h1], {k: e[k][h1] for k in keep}, T)['chi_total']
            r['chi_half2'] = fdt.susceptibility(obs[m][h2], {k: e[k][h2] for k in keep}, T)['chi_total']
            rows.append({**base, 'metric': m, **r})
        f32 = lambda x: np.asarray(x, dtype=np.float32)  # noqa: E731
        np.savez_compressed(os.path.join(npz_dir, f'{uid}.npz'), T=T, sequence=seq, labels=np.array(labels),
                            frame_index=info['production_start'] + info['stride'] * np.arange(n),
                            E_per_force=f32(E0), dEdT_per_force=f32(dEdT),
                            **{k: f32(v) for k, v in e.items()}, **{m: f32(a) for m, a in obs.items()})
        return uid, rows, None, round(time.time() - t0, 1)
    except Exception as ex:
        return uid, None, f'{ex!r} | {traceback.format_exc()[-400:]}'.replace('\n', ' '), round(time.time() - t0, 1)


def run_dirs(root: str, T: float, ids_file: str = None) -> dict:
    tdir = os.path.join(root, f'{int(T)}K')
    names = [l.strip() for l in open(ids_file) if l.strip()] if ids_file else sorted(os.listdir(tdir))
    return {u: os.path.join(tdir, u) for u in names
            if os.path.isfile(os.path.join(tdir, u, f'{u}.dcd')) and os.path.isfile(os.path.join(tdir, u, 'top.pdb'))}


def check(root: str, T: float, name: str = None, ionic_M: float = IONIC_M, platform: str = 'CPU', n_frames: int = 20):
    import openmm as mm
    dirs = run_dirs(root, T)
    name = name or next(iter(dirs))
    xyz, seq, info = load_run(dirs[name], name)
    xyz = xyz[:n_frames]
    print(f'{name}: N {len(seq)}, DCD frames {info["n_frames_total"]} (expected {bc.expected_frames(len(seq))}), '
          f'max bond {info["max_bond_raw_nm"]:.3f} nm raw / {info["max_bond_nm"]:.3f} nm used\n')
    s0, s1 = build_system(name, seq, T, ionic_M, xyz[0]), build_system(name, seq, T + DT, ionic_M, xyz[0])
    print('forces:')
    for d in describe_forces(s0):
        print(f"  [{d['index']}] {d['class']:24s} label={d['label']:7s} cutoff={d['cutoff_nm']} globals={d['globals']}")
        if d['expression']:
            print(f"       {d['expression']}")
    a, b = mm.XmlSerializer.serialize(s0).splitlines(), mm.XmlSerializer.serialize(s1).splitlines()
    diff = [l for l in difflib.unified_diff(a, b, lineterm='', n=0) if l[:1] in '+-' and l[:3] not in ('+++', '---')]
    print(f'\nserialized system, T vs T + {DT} K: {len(diff) // 2} changed lines; first ones:')
    for l in diff[:12]:
        print('   ', l[:200])
    T_ = float(T)
    expect = -0.5 * (1 / T_ + fdt.deps_w_dT(T_) / fdt.eps_w(T_)) * DT
    print(f'paper formula (kappa from salt and eps_r(T)): d ln kappa = {expect:+.6f} for +{DT} K')
    for d0, d1 in zip(describe_forces(s0), describe_forces(s1)):
        for g, v in d0['globals'].items():
            if d1['globals'].get(g) != v and v:
                print(f"  global '{g}' in force [{d0['index']}]: {v:.6g} -> {d1['globals'][g]:.6g} "
                      f"(d ln = {np.log(d1['globals'][g] / v):+.6f})")
    t0 = time.time()
    labels, E0, dEdT = energies_at(name, seq, xyz, T, ionic_M, platform)
    t_e = time.time() - t0
    _, _, dEdT_half = energies_at(name, seq, xyz, T, ionic_M, platform, dT=DT / 2)
    ev = Evaluator(build_system(name, seq, T, ionic_M, xyz[0]), platform)
    tot = np.array([ev.total(x) for x in xyz])
    print(f'\ngroup energies vs total: max |sum - total| = {np.abs(E0.sum(1) - tot).max():.2e} kJ/mol')
    print('mean energy (kJ/mol) and mean dE/dT (kJ/mol/K) per force:')
    for k, lab in enumerate(labels):
        print(f'  [{k}] {lab:7s} E {E0[:, k].mean():12.4f}   dE/dT {dEdT[:, k].mean():+.6e}   '
              f'(DT/2: {dEdT_half[:, k].mean():+.6e})')
    print(f'\ntiming: {t_e / len(xyz) * K_FRAMES / 60:.1f} min per protein for the energies of {K_FRAMES} frames (3 systems)')


def run(root, T, out, ids_file=None, workers=4, platform='CPU', limit=None, ionic_M=IONIC_M):
    os.makedirs(out, exist_ok=True)
    npz_dir = os.path.join(out, f'npz_{int(T)}K')
    os.makedirs(npz_dir, exist_ok=True)
    os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    csv_path = os.path.join(out, f'recharged_fdt_{int(T)}K.csv')
    fail_path = os.path.join(out, 'logs', f'failed_{int(T)}K.tsv')
    dirs = run_dirs(root, T, ids_file)
    done = set(pd.read_csv(csv_path, usecols=['UniProt_ID'])['UniProt_ID']) if os.path.exists(csv_path) else set()
    todo = [u for u in dirs if u not in done][:limit]
    print(f'{T:g} K: {len(dirs)} runs, {len(done)} done, {len(todo)} to run, {workers} workers on {platform}', flush=True)
    t0, n_fail = time.time(), 0
    with mp.get_context('spawn').Pool(workers) as pool:
        tasks = [(u, dirs[u], float(T), ionic_M, platform, npz_dir) for u in todo]
        for i, (uid, rows, err, sec) in enumerate(pool.imap_unordered(analyse, tasks), 1):
            if rows:
                pd.DataFrame(rows).to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)
            else:
                n_fail += 1
                with open(fail_path, 'a') as fh:
                    fh.write(f'{uid}\t{err}\n')
            if i % 25 == 0 or err:
                el = (time.time() - t0) / 3600
                print(f'[{T:g} K {i}/{len(todo)}] {el:.2f} h, ~{el / i * (len(todo) - i):.1f} h left, {n_fail} failed'
                      + (f' | {uid}: {err[:160]}' if err else ''), flush=True)
    bc.write_provenance(csv_path, __file__, {'temperature_K': T, 'ionic_M': ionic_M, 'dT_central_difference_K': DT,
                                             'frames': f'{K_FRAMES} evenly spaced production frames', 'equilibration_fraction': EQUIL_FRAC,
                                             'energies': 'OpenMpipi get_mpipi_system at T and T +/- DT, one force group per force',
                                             'platform': platform, 'workers': workers}, [root])
    print(f'{T:g} K finished: {csv_path}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    c = sub.add_parser('check')
    c.add_argument('--runs', required=True, help='root holding <T>K/<ID>/ (run_openmpipi_idp.py --outdir)')
    c.add_argument('--temp', required=True, type=float)
    c.add_argument('--name', default=None)
    c.add_argument('--platform', default='CPU', choices=['CPU', 'CUDA', 'Reference'])
    c.add_argument('--ionic', type=float, default=IONIC_M)
    lc = sub.add_parser('logcheck')
    lc.add_argument('--runs', required=True)
    lc.add_argument('--temp', required=True, type=float)
    lc.add_argument('--name', default=None)
    lc.add_argument('--platform', default='CPU', choices=['CPU', 'CUDA', 'Reference'])
    lc.add_argument('--ionic', type=float, default=IONIC_M)
    r = sub.add_parser('run')
    r.add_argument('--runs', required=True)
    r.add_argument('--temp', required=True, type=float)
    r.add_argument('--out', required=True)
    r.add_argument('--ids', default=None)
    r.add_argument('--workers', type=int, default=4)
    r.add_argument('--platform', default='CPU', choices=['CPU', 'CUDA', 'Reference'])
    r.add_argument('--limit', type=int, default=None)
    r.add_argument('--ionic', type=float, default=IONIC_M)
    a = ap.parse_args()
    if a.mode == 'check':
        check(a.runs, a.temp, a.name, a.ionic, a.platform)
    elif a.mode == 'logcheck':
        logcheck(a.runs, a.temp, a.name, a.ionic, a.platform)
    else:
        run(a.runs, a.temp, a.out, a.ids, a.workers, a.platform, a.limit, a.ionic)
