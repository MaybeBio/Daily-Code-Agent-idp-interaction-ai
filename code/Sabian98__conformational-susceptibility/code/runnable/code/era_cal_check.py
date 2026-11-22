"""
era_cal_check.py — check of the GAIVI_ERACAL inputs (make_era_calibrated_inputs.py) before any run: on frames of the
plain 277 K ER-NTD runs (GAIVI_ERA), the energy of the calibrated system minus that of the plain system must equal the
bias energy of the reweighting, sum_c gamma_c E_c, frame by frame, for every form and both variants.
  CALVADOS-2       a CALVADOS system built at 277 K and the run's ionic strength (config and components of a held-out run,
                   with the form's FASTA) with the GAIVI_ERA residues file and with the package's residues_<form>_<v>_277K.csv;
                   reference E_c from tea_reweight.class_energies and tea_calibrate.calvados_hyd_energies on the sequence
                   with the phosphosite read as Ser (era_direct.py)
  Mpipi-Recharged  mpipi_recharged_fdt.build_system (+ era_direct.apply_phospho for pS118), without and with
                   run_openmpipi_cal.add_hydration_forcing (the package copy) and the package's delta file; reference E_c
                   from tea_calibrate.recharged_class_energies on the run sequence (Glu at the phosphosite)
S1cap (capped, many-body): the added CustomGBForce against saturating_forcing's capped energy from per-residue contacts.
60 frames of the second half of each plain 277 K run (0.195 M runs for WT, S118D, pS118; 0.095 M for WT, pS118).
    BENDER_REANALYSIS_DIR=... python era_cal_check.py [--pkg $CS_DATA/Scripts/GAIVI_ERACAL]    (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/era/era_cal_check.txt
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
import era_direct as ED  # noqa: E402
import heldout_analysis as HA  # noqa: E402
import heldout_exploratory as HX  # noqa: E402
import make_era_calibrated_inputs as ME  # noqa: E402
import make_pnt_calibrated_inputs as MP  # noqa: E402
import saturating_forcing as SF  # noqa: E402

TCal, W = HA.TCal, HA.W
CFG_SRC = f'{HA.BM}/heldout_CALVADOS2/303K/HO-PNt-I175-pH74'
NF = 60
OUTF = f'{ED.OUTROOT}/era_cal_check.txt'


def bias_gam(model, variant):
    gt, gh = MP.gam(model, variant)
    return {'tea': gt, 'hyd_all': gh} if variant == 'A0' else HX.chg_gam(gh, gt)


def frames(model, form, ionic):
    tag = f'I{round(ionic * 1000)}'
    name = f'ERa-{form}-{tag}-pH74' if model == 'calvados' else f'ERa-{form}-{tag}'
    d = os.path.join(ED.CFG[model]['root'], '277K', name)
    t = md.load(os.path.join(d, f'{name}.dcd'), top=os.path.join(d, 'top.pdb'))
    n = t.n_frames
    x = t.xyz[n // 2 + np.linspace(0, n - n // 2 - 1, NF).astype(int)].astype(float)
    if HA.bc.max_bond(x) > 1.0:
        x = HA.bc.unwrap_chain(x, HA.box_nm(187))
    return x


def calvados_energies(seq, ionic, residues, xyz, add=None):
    import openmm as mm
    from calvados.sim import Sim
    cfg = yaml.safe_load(open(os.path.join(CFG_SRC, 'config.yaml')))
    comp = yaml.safe_load(open(os.path.join(CFG_SRC, 'components.yaml')))
    cfg.update(temp=float(ME.T), ionic=float(ionic), pH=ME.PH, platform='CPU', verbose=False, sysname='ER')
    cfg.pop('runtime', None)
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, 'ER.fasta')
        open(fa, 'w').write(f'>ER\n{seq}\n')
        comp['defaults'].update(ffasta=fa, fresidues=residues)
        comp['system'] = {'ER': {}}
        sim = Sim(tmp, cfg, comp)
        sim.build_system()
        if add is not None:
            add(sim.system)
        groups_ok(sim.system)
        ctx = mm.Context(sim.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        out = []
        for x in xyz:
            ctx.setPositions(x)
            out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def recharged_energies(seq, ionic, xyz, phospho, delta=None, wrapper=None):
    import openmm as mm
    import mpipi_recharged_fdt as R
    system = R.build_system('ER', seq, float(ME.T), ionic, xyz[0])
    if phospho != '-':
        ED.apply_phospho(system, phospho, seq)
    if delta is not None:
        kind, forcing = wrapper.read_forcing(delta)
        wrapper.add_hydration_forcing(system, forcing) if kind == 'linear' else wrapper.add_capped_forcing(system, *forcing)
    groups_ok(system)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
    out = []
    for x in xyz:
        ctx.setPositions(x)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def groups_ok(system):
    """OpenMM on GPUs requires one cutoff per force group (the first PNt runs failed on this)."""
    import openmm as mm
    g = {}
    for f in system.getForces():
        if isinstance(f, (mm.CustomNonbondedForce, mm.CustomGBForce, mm.NonbondedForce)):
            c = f.getCutoffDistance()
            c = c.value_in_unit(mm.unit.nanometer) if hasattr(c, 'value_in_unit') else float(c)
            g.setdefault(f.getForceGroup(), set()).add(round(c, 9))
    bad = {k: v for k, v in g.items() if len(v) > 1}
    if bad:
        raise RuntimeError(f'force groups with more than one cutoff: {bad}')


def capped_reference(model, xyz, seq, T):
    """gamma_TEA E_tsat + gamma_h E_hsat of saturating_forcing.py (charged-group form, c_sat 2)."""
    gt, gh = ME.capped_gamma(model)
    c = SF.contacts(model, xyz, seq)
    wt, wh = SF.weights(model, seq, float(T), 'DEKR')
    return gt * SF.energy(model, c, wt, ME.CSAT) + gh * SF.energy(model, c, wh, ME.CSAT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pkg', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_ERACAL')))
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location('run_openmpipi_cal', os.path.join(a.pkg, 'mpipi', 'run_openmpipi_cal.py'))
    wrapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wrapper)
    cspec = importlib.util.spec_from_file_location('run_calvados_cal', os.path.join(a.pkg, 'calvados', 'run_calvados_cal.py'))
    cwrap = importlib.util.module_from_spec(cspec)
    cspec.loader.exec_module(cwrap)
    cfa = ME.read_fasta(os.path.join(a.pkg, 'calvados', 'forms.fasta'))
    mfa = ME.read_fasta(os.path.join(a.pkg, 'mpipi', 'forms.fasta'))
    L = ['ER-NTD at 277 K: direct calibrated energy minus plain energy against the reweighting bias (kJ/mol), '
         f'{NF} plain frames per row', '']
    worst = 0.0
    for model in ('calvados', 'recharged'):
        for form, ionic in ME.RUNS:
            xyz = frames(model, form, ionic)
            phos = f'{ME.SITE}:2' if form == 'pS118' else '-'
            ref_seq = None
            if model == 'calvados':
                seq = cfa[form]
                ref_seq = seq.replace('B', 'S')
                E = W.class_energies(xyz, ref_seq, float(ME.T))
                E.update(TCal.calvados_hyd_energies(xyz, ref_seq, float(ME.T)))
                E0 = calvados_energies(seq, ionic, f'{ME.ERA}/calvados/residues_{form}.csv', xyz)
            else:
                seq = mfa[form]
                E = TCal.recharged_class_energies(xyz, seq, float(ME.T), HA.box_nm(len(seq)))
                E0 = recharged_energies(seq, ionic, xyz, phos)
            z = {f'E_{k}': v for k, v in E.items()}
            for variant in ('A0', 'S1', 'S1cap'):
                if variant == 'S1cap':
                    b = capped_reference(model, xyz, ref_seq if model == 'calvados' else seq, ME.T)
                else:
                    b = HA.bias(z, bias_gam(model, variant), float(ME.T)) * HA.R_GAS * ME.T
                if model == 'calvados' and variant == 'S1cap':
                    w, csat = cwrap.read_cap(os.path.join(a.pkg, 'calvados', f'w_S1cap_{ME.T}K.tsv'))
                    E1 = calvados_energies(seq, ionic, os.path.join(a.pkg, 'calvados', f'residues_{form}_plain.csv'), xyz,
                                           add=lambda sysm: cwrap.add_capped_forcing_calvados(sysm, w, csat, seq))
                elif model == 'calvados':
                    E1 = calvados_energies(seq, ionic, os.path.join(a.pkg, 'calvados', f'residues_{form}_{variant}_{ME.T}K.csv'), xyz)
                else:
                    fname = f'w_S1cap_{ME.T}K.tsv' if variant == 'S1cap' else f'delta_{variant}_{ME.T}K.tsv'
                    E1 = recharged_energies(seq, ionic, xyz, phos, os.path.join(a.pkg, 'mpipi', fname), wrapper)
                dE = E1 - E0
                err = np.abs(dE - b).max()
                rel = err / max(np.abs(b).max(), 1e-12)
                worst = max(worst, rel)
                L.append(f'{HA.LABEL[model]:16s} {form:6s} {ionic} M {variant}: bias mean {b.mean():+8.3f} sd {b.std():6.3f}; '
                         f'max |direct - bias| {err:.2e} ({rel:.1e} of max |bias|)')
                print(L[-1], flush=True)
    L += ['', f'largest relative difference {worst:.1e}: ' + ('PASS (< 1e-4)' if worst < 1e-4 else 'FAIL')]
    open(OUTF, 'w').write('\n'.join(L) + '\n')
    print(L[-1])


if __name__ == '__main__':
    main()
