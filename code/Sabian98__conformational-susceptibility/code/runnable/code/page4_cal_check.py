"""
page4_cal_check.py — check of the GAIVI_PAGE4CAL inputs (make_page4_calibrated_inputs.py) before any run: on frames of
the plain PAGE4 heating runs (page4_direct.py) at each temperature, the energy of the calibrated system minus that of the
plain system must equal the saturating hydration energy of the analysis, gamma_TEA E_tsat + gamma_h E_hsat
(saturating_forcing.contacts / weights / energy, charged-group set, c_sat 2), frame by frame, for every form.
  CALVADOS-2       a CALVADOS system built from the plain run's config.yaml and components.yaml at its temperature, with the
                   package's plain residues file, without and with run_calvados_cal.add_capped_forcing_calvados (package
                   copy, w_cal_<T>K.tsv); reference on the sequence with each phosphoresidue read as its Ser or Thr
  Mpipi-Recharged  mpipi_recharged_fdt.build_system (+ era_direct.apply_phospho, the verbatim copy of
                   run_openmpipi_page4.apply_phospho) without and with run_openmpipi_cal.add_capped_forcing (package copy);
                   reference on the run sequence (Glu beads at phosphosites)
Every force group must have one cutoff (GPU requirement). 30 frames of the second half of each plain run.
    BENDER_REANALYSIS_DIR=... python page4_cal_check.py [--pkg $CS_DATA/Scripts/GAIVI_PAGE4CAL]    (openmpipi env)
Output: heldout_temperature/capped_forcing/page4/page4_cal_check.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import os  # noqa: E402
import tempfile  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import era_cal_check as EC  # noqa: E402
import era_direct as ED  # noqa: E402
import make_era_calibrated_inputs as ME  # noqa: E402
import make_page4_calibrated_inputs as MPC  # noqa: E402
import page4_direct as PD  # noqa: E402
import saturating_forcing as SF  # noqa: E402

NF = 30
OUTF = os.path.join(PD.BM, 'heldout_temperature', 'capped_forcing', 'page4', 'page4_cal_check.txt')


def load(module, path):
    spec = importlib.util.spec_from_file_location(module, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def frames(d, name, natoms_box=None):
    t = md.load(os.path.join(d, f'{name}.dcd'), top=os.path.join(d, 'top.pdb'))
    n = t.n_frames
    x = t.xyz[n // 2 + np.linspace(0, n - n // 2 - 1, NF).astype(int)].astype(float)
    if PD.HC.bc.max_bond(x) > 1.0:
        x = PD.HC.bc.unwrap_chain(x, natoms_box or 50.0)
    return x


def calvados_energies(d, name, T, residues, xyz, add=None):
    import openmm as mm
    from calvados.sim import Sim
    cfg = yaml.safe_load(open(os.path.join(d, 'config.yaml')))
    comp = yaml.safe_load(open(os.path.join(d, 'components.yaml')))
    cfg.update(temp=float(T), platform='CPU', verbose=False)
    cfg.pop('runtime', None)
    comp['defaults'].update(ffasta=os.path.join(d, f'{name}.fasta'), fresidues=residues)
    with tempfile.TemporaryDirectory() as tmp:
        sim = Sim(tmp, cfg, comp)
        sim.build_system()
        if add is not None:
            add(sim.system)
        EC.groups_ok(sim.system)
        ctx = mm.Context(sim.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        out = []
        for x in xyz:
            ctx.setPositions(x)
            out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def recharged_energies(seq, T, ionic, spec, xyz, wfile=None, wrapper=None):
    import openmm as mm
    import mpipi_recharged_fdt as R
    system = R.build_system('PAGE4', seq, float(T), ionic, xyz[0])
    if spec not in ('-', '', None):
        ED.apply_phospho(system, spec, seq)
    if wfile is not None:
        kind, forcing = wrapper.read_forcing(wfile)
        assert kind != 'linear', 'the PAGE4CAL weight files must be the saturating form'
        wrapper.add_capped_forcing(system, *forcing)
    EC.groups_ok(system)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
    out = []
    for x in xyz:
        ctx.setPositions(x)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pkg', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4CAL')))
    a = ap.parse_args()
    cwrap = load('run_calvados_cal', os.path.join(a.pkg, 'calvados', 'run_calvados_cal.py'))
    mwrap = load('run_openmpipi_cal', os.path.join(a.pkg, 'mpipi', 'run_openmpipi_cal.py'))
    L = ['PAGE4 heating series: direct calibrated energy minus plain energy against the saturating energy of the analysis '
         f'(kJ/mol), {NF} plain frames per row', '']
    worst = 0.0
    for model, sub in (('calvados', 'calvados'), ('recharged', 'mpipi')):
        runs = pd.read_csv(os.path.join(a.pkg, sub, 'runs.tsv'), sep='\t')
        fa = ME.read_fasta(os.path.join(a.pkg, sub, 'forms.fasta'))
        letters = {}                                     # per form (letters are assigned form by form)
        if model == 'calvados':
            for f in fa:
                letters[f] = MPC.phospho_base(os.path.join(a.pkg, sub, f'residues_{f}_plain.csv'))
        for _, r in runs.iterrows():
            T, f, seq = int(r.temp), r.form, fa[r.form]
            if model == 'calvados':
                name = f'PAGE4-{f}-I150-pH75'
                d = os.path.join(PD.MODELS[model]['root'], f'{T}K', name)
                run_seq = ''.join(ME.read_fasta(os.path.join(d, f'{name}.fasta')).values())
                assert run_seq == seq, f'{name}: package sequence differs from the run'
                xyz = frames(d, name)
                res = os.path.join(a.pkg, sub, r.residues)
                w, csat = cwrap.read_cap(os.path.join(a.pkg, sub, r.cap))
                E0 = calvados_energies(d, name, T, res, xyz)
                E1 = calvados_energies(d, name, T, res, xyz, add=lambda s: cwrap.add_capped_forcing_calvados(s, w, csat, seq))
                ref_seq = ''.join(letters[f].get(c, c) for c in seq)
            else:
                import mpipi_recharged_fdt as R
                name = f'PAGE4-{f}-I150'
                d = os.path.join(PD.MODELS[model]['root'], f'{T}K', name)
                xyz = frames(d, name, R.box_nm(len(seq)))
                top_seq = ''.join(x.code for x in md.load_topology(os.path.join(d, 'top.pdb')).residues)
                assert top_seq == seq, f'{name}: package sequence differs from the run'
                E0 = recharged_energies(seq, T, float(r.ionic), r.phospho, xyz)
                E1 = recharged_energies(seq, T, float(r.ionic), r.phospho, xyz, os.path.join(a.pkg, sub, r.delta), mwrap)
                ref_seq = seq
            b = EC.capped_reference(model, xyz, ref_seq, T)
            dE = E1 - E0
            err = np.abs(dE - b).max()
            rel = err / max(np.abs(b).max(), 1e-12)
            worst = max(worst, rel)
            L.append(f'{PD.MODELS[model]["label"]:16s} {f:9s} {T} K: saturating energy mean {b.mean():+8.3f} sd {b.std():6.3f}; '
                     f'max |direct - analysis| {err:.2e} ({rel:.1e} of max)')
            print(L[-1], flush=True)
    L += ['', f'largest relative difference {worst:.1e}: ' + ('PASS (< 1e-4)' if worst < 1e-4 else 'FAIL')]
    os.makedirs(os.path.dirname(OUTF), exist_ok=True)
    open(OUTF, 'w').write('\n'.join(L) + '\n')
    print(L[-1])


if __name__ == '__main__':
    main()
