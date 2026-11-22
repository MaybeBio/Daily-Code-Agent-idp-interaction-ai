"""
pnt_cal_check.py — check of the GAIVI_PNTCAL inputs (make_pnt_calibrated_inputs.py) before any run: on frames of the
held-out PNt runs, the energy of the calibrated system minus that of the plain system must equal the bias energy of the
reweighting, sum_c gamma_c E_c (heldout_analysis.bias x RT, from the stored held-out npz), frame by frame.
  CALVADOS-2       the run's own CALVADOS system (config.yaml, components.yaml of the held-out run) built with the
                   standard residues file and with the package's residues_<variant>_<T>K.csv
  Mpipi-Recharged  mpipi_recharged_fdt.build_system, without and with run_openmpipi_cal.add_hydration_forcing (the
                   package copy) and the package's delta_<variant>_<T>K.tsv
Sets and temperatures checked: PNt-I175 at 283 and 331 K and PNt-I025 at 331 K, A0 and S1, 60 stored frames each.
    BENDER_REANALYSIS_DIR=... python pnt_cal_check.py [--pkg $CS_DATA/Scripts/GAIVI_PNTCAL]    (openmpipi env)
Output: heldout_temperature/results/exploratory/pnt_cal_check.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import os  # noqa: E402
import tempfile  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import heldout_analysis as HA  # noqa: E402
import heldout_exploratory as HX  # noqa: E402
import make_heldout_inputs as MH  # noqa: E402
import make_pnt_calibrated_inputs as MP  # noqa: E402

CHECKS = [('PNt-I175', 283), ('PNt-I175', 331), ('PNt-I025', 331)]
NF = 60
RUNDIR = {'calvados': lambda st, T: f'{HA.BM}/heldout_CALVADOS2/{T}K/HO-{st}-pH74',
          'recharged': lambda st, T: f'{HA.BM}/heldout_MPIPI_Recharged/{T}K/HO-{st}'}
NPZ = {'calvados': lambda st, T: f'{HA.OUT}/npz_calvados/HO-{st}-pH74_{T}K.npz',
       'recharged': lambda st, T: f'{HA.OUT}/npz_recharged/HO-{st}_{T}K.npz'}


def bias_gam(model, variant):
    gt, gh = MP.gam(model, variant)
    return {'tea': gt, 'hyd_all': gh} if variant == 'A0' else HX.chg_gam(gh, gt)


def frames(model, st, T):
    z = dict(np.load(NPZ[model](st, T)))
    d = RUNDIR[model](st, T)
    name = os.path.basename(d)
    t = md.load(os.path.join(d, f'{name}.dcd'), top=os.path.join(d, 'top.pdb'))
    if t.n_frames != int(z['n_frames_total']):
        raise ValueError(f'{name} {T} K: DCD has {t.n_frames} frames, npz lists {z["n_frames_total"]}')
    pick = np.linspace(0, len(z['frame_index']) - 1, NF).astype(int)
    xyz = t.xyz[z['frame_index'][pick]].astype(float)
    ref = {k: v[pick] for k, v in z.items() if k.startswith('E_')}
    return d, name, xyz, ref


def calvados_energies(d, name, T, residues, xyz):
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
        ctx = mm.Context(sim.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        out = []
        for x in xyz:
            ctx.setPositions(x)
            out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def recharged_energies(seq, T, ionic, xyz, delta=None, wrapper=None):
    import openmm as mm
    import mpipi_recharged_fdt as R
    system = R.build_system('PNt', seq, T, ionic, xyz[0])
    if delta is not None:
        wrapper.add_hydration_forcing(system, delta)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
    out = []
    for x in xyz:
        ctx.setPositions(x)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pkg', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PNTCAL')))
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location('run_openmpipi_cal', os.path.join(a.pkg, 'mpipi', 'run_openmpipi_cal.py'))
    wrapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wrapper)
    L = ['Direct calibrated energy minus plain energy against the reweighting bias sum_c gamma_c E_c (kJ/mol), '
         f'{NF} stored held-out frames per row', '']
    worst = 0.0
    for model in ('calvados', 'recharged'):
        for st, T in CHECKS:
            d, name, xyz, ref = frames(model, st, T)
            ionic = MP.SETS[st][0]
            if model == 'calvados':
                E0 = calvados_energies(d, name, T, MH.RESIDUES, xyz)
            else:
                E0 = recharged_energies(MH.PNT, float(T), ionic, xyz)
            for variant in ('A0', 'S1'):
                b = HA.bias(ref, bias_gam(model, variant), float(T)) * HA.R_GAS * T
                if model == 'calvados':
                    E1 = calvados_energies(d, name, T, os.path.join(a.pkg, 'calvados', f'residues_{variant}_{T}K.csv'), xyz)
                else:
                    delta = wrapper.read_delta(os.path.join(a.pkg, 'mpipi', f'delta_{variant}_{T}K.tsv'))
                    E1 = recharged_energies(MH.PNT, float(T), ionic, xyz, delta, wrapper)
                dE = E1 - E0
                err = np.abs(dE - b).max()
                rel = err / max(np.abs(b).max(), 1e-12)
                worst = max(worst, rel)
                L.append(f'{HA.LABEL[model]:16s} {st} {T} K {variant}: bias mean {b.mean():+9.3f} sd {b.std():7.3f}; '
                         f'direct - bias max |diff| {err:.2e} kJ/mol ({rel:.1e} of max |bias|); r {np.corrcoef(dE, b)[0, 1]:.6f}')
                print(L[-1], flush=True)
    L += ['', f'largest relative difference {worst:.1e}: ' + ('PASS (< 1e-4)' if worst < 1e-4 else 'FAIL')]
    os.makedirs(os.path.join(HA.OUT, 'exploratory'), exist_ok=True)
    open(os.path.join(HA.OUT, 'exploratory', 'pnt_cal_check.txt'), 'w').write('\n'.join(L) + '\n')
    print(L[-1])


if __name__ == '__main__':
    main()
