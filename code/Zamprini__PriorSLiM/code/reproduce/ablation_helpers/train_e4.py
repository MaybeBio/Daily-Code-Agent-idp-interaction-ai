from __future__ import annotations
import argparse
import gc
import logging
import os
from pathlib import Path
import sys
import time
os.environ.setdefault('PYTORCH_ENABLE_MPS_FALLBACK', '1')
from e4_data import HERE, ExperimentData, csr_events, read, require
from train_jwarm import batches, device_for, digest, identity, load_pt, metrics, save_pt, seed_all, write_json, jsonable
from train_jwarm import cpu_tree, exclusive, fit_stage, release, rng_state, restore_rng
from model_e4 import BinaryGate, DecoderReadout, profile_sum
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
import torch
from torch.nn import functional as F

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['check', 'decoder', 'profile', 'controls', 'refit'])
    p.add_argument('--project-root', type=Path, default=HERE.parents[1])
    p.add_argument('--e1-results', type=Path, default=HERE / 'E1/results')
    p.add_argument('--proximal-results', type=Path)
    p.add_argument('--output-root', type=Path, default=HERE / 'E4/results')
    p.add_argument('--store-dir', type=Path)
    p.add_argument('--decoder-checkpoint', type=Path)
    p.add_argument('--variant', choices=['A', 'B', 'C'], default='C')
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--device', choices=['cpu', 'mps', 'cuda', 'auto'], default='auto')
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--readout-epochs', type=int, default=20)
    p.add_argument('--readout-patience', type=int, default=4)
    p.add_argument('--readout-lr', type=float, default=0.0003)
    p.add_argument('--readout-wd', type=float, default=0.0001)
    p.add_argument('--gate-epochs', type=int, default=100)
    p.add_argument('--gate-patience', type=int, default=15)
    p.add_argument('--gate-lr', type=float, default=0.03)
    p.add_argument('--gate-lambdas', nargs='+', type=float, default=[0.0, 0.03, 0.1, 0.3, 1.0, 3.0])
    p.add_argument('--budgets', nargs='+', type=int, default=[64, 128, 256, 512, 1024, 2048, 4096, 5845])
    p.add_argument('--random-repeats', type=int, default=10)
    p.add_argument('--decoder-chunk-nnz', type=int, default=2048)
    p.add_argument('--refit-ridge', type=float, default=0.0001)
    p.add_argument('--refit-max-iter', type=int, default=3000)
    p.add_argument('--refit-tol', type=float, default=1e-05)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args(argv)
    for k in ['project_root', 'e1_results', 'output_root']:
        setattr(args, k, getattr(args, k).expanduser().resolve())
    args.proximal_results = (args.proximal_results or args.e1_results).resolve()
    args.store_dir = (args.store_dir or args.project_root / 'feature_cache/motifae_sparse_windows_idr_coverage').resolve()
    default_checkpoint = Path(os.environ.get(
        'MOTIFAE_CHECKPOINT', args.project_root / 'motifae/MotifAE_k3_step_80000.pt'
    ))
    args.decoder_checkpoint = (args.decoder_checkpoint or default_checkpoint).expanduser().resolve()
    args.gate_lambdas = sorted(set(args.gate_lambdas))
    args.budgets = sorted(set(args.budgets))
    if any((not np.isfinite(x) or x < 0 for x in args.gate_lambdas)) or not args.gate_lambdas:
        p.error('gate lambdas must be finite and nonnegative')
    if any((k < 1 or k > 5845 for k in args.budgets)):
        p.error('budgets must be within 1..5845')
    for key in ['batch_size', 'threads', 'readout_epochs', 'readout_patience', 'gate_epochs', 'gate_patience', 'random_repeats', 'decoder_chunk_nnz', 'refit_max_iter']:
        if getattr(args, key) <= 0:
            p.error(f'{key} must be positive')
    for key in ['readout_lr', 'gate_lr', 'refit_ridge', 'refit_tol']:
        if not np.isfinite(getattr(args, key)) or getattr(args, key) <= 0:
            p.error(f'{key} must be finite and positive')
    if not np.isfinite(args.readout_wd) or args.readout_wd < 0:
        p.error('readout_wd must be finite and nonnegative')
    return args

def training_contract(args):
    keys = ['variant', 'seed', 'batch_size', 'readout_epochs', 'readout_patience', 'readout_lr', 'readout_wd', 'gate_epochs', 'gate_patience', 'gate_lr', 'gate_lambdas', 'budgets', 'random_repeats', 'decoder_chunk_nnz', 'refit_ridge', 'refit_max_iter', 'refit_tol']
    return {k: getattr(args, k) for k in keys}

def record(out, data, point_id, method, family, coefficient, logits, **extra):
    out.mkdir(parents=True, exist_ok=True)
    coefficient = np.asarray(coefficient, float)
    support = np.flatnonzero(coefficient > 1e-10)
    result = dict(point_id=point_id, method=method, family=family, variant=data.args.variant, seed=data.args.seed, K=len(support), candidate_count=len(data.ids), feature_fraction=len(support) / len(data.ids), validation_only=True, test_evaluated=False, eligible=True)
    result.update(extra)
    for split, z in logits.items():
        m = metrics(data.splits[split].labels, z, data.splits[split].negative_types)
        result.update({f'{split}_{k}': v for k, v in m.items()})
        if split == 'val':
            (out / 'predictions').mkdir(exist_ok=True)
            np.savez_compressed(out / f'predictions/{point_id}.npz', logits=z, labels=data.splits[split].labels, window_ids=data.splits[split].window_ids.astype(str))
    (out / 'solutions').mkdir(exist_ok=True)
    save_pt(out / f'solutions/{point_id}.pt', dict(feature_ids=torch.from_numpy(data.ids), coefficient=torch.from_numpy(coefficient), selected_feature_ids=torch.from_numpy(data.ids[support]), meaning='gate multiplier on original profile evidence' if family == 'profile' else 'binary latent gate'))
    write_json(out / f'points/{point_id}.json', result)
    logging.info('%s K=%d valAP=%.6f', point_id, len(support), result['val_ap'])
    return result

def profile_logits(data, multipliers):
    return {s: data.explicit[s] + data.evidence[s] @ multipliers for s in data.splits}

def predict_head(data, head, args, device, gates=None, fused=False):
    result = {}
    head.eval()
    with torch.inference_mode():
        for s, split in data.splits.items():
            z = np.empty(len(split), float)
            for take in batches(len(split), args.batch_size):
                r = head(data.raw.events(s, take, device), gates, fused=fused)
                z[take] = data.explicit[s][take] + r.cpu().numpy()
            result[s] = z
    return result

def profile_reference(args, data, out):
    if not (out / 'points/explicit_only.json').exists():
        record(out, data, 'explicit_only', 'explicit_only', 'profile', np.zeros(len(data.ids)), data.explicit, reference_point='profile_full')
    point = out / 'points/profile_full.json'
    if not point.exists():
        record(out, data, 'profile_full', 'full', 'profile', np.ones(len(data.ids)), profile_logits(data, np.ones(len(data.ids))), reference_point='profile_full', positional_parameters=len(data.ids) * 29)

def decoder_stage(args, data, out, signature, device):
    seed_all(args.seed + 400)
    head = DecoderReadout(data.vectors, data.bias, args.decoder_chunk_nnz).to(device)
    forward = lambda model, split, take: torch.tensor(data.explicit[split][take], dtype=torch.float32, device=device) + model(data.raw.events(split, take, device))
    checkpoint = fit_stage(head, forward, data.splits, out / 'decoder_readout', args, device, stage='E4 decoder all-open readout', epochs=args.readout_epochs, patience=args.readout_patience, lr=args.readout_lr, wd=args.readout_wd, stage_signature=identity([signature, 'decoder_readout']))
    count = sum((p.numel() for p in head.parameters()))
    head.freeze_for_gates()
    take = np.arange(min(4, len(data.splits['val'])))
    batch = data.raw.events('val', take, device)
    with torch.no_grad():
        plain = head(batch)
        fast = head(batch, fused=True)
        require(torch.allclose(plain, fast, atol=2e-05, rtol=2e-05), 'Frozen decoder+Linear composition mismatch')
    write_json(out / 'decoder_composition_check.json', dict(max_error=float((plain - fast).abs().max().cpu()), raw_decoder_shape=[len(data.ids), 1280], bias_included=True, trainable_decoder_parameters=0, readout_parameters=count, frozen_after_readout=True))
    if not (out / 'points/decoder_full.json').exists():
        record(out, data, 'decoder_full', 'full', 'decoder', np.ones(len(data.ids)), predict_head(data, head, args, device, fused=True), reference_point='decoder_full', best_epoch=checkpoint['best_epoch'], readout_parameters=count, decoder_weight_parameters=len(data.ids) * 1280, decoder_bias_parameters=1280)
        record(out, data, 'decoder_zero', 'zero_gate', 'decoder', np.zeros(len(data.ids)), predict_head(data, head, args, device, gates=torch.zeros(len(data.ids), device=device), fused=True), reference_point='decoder_full', readout_parameters=count, note='H(decoder bias) may be nonzero; includes padding mask effects, not identical to explicit-only')
    for lam in args.gate_lambdas:
        fit_gate(args, data, out, signature, device, 'decoder', lam, head)
    del head
    release(device)

def fit_gate(args, data, out, signature, device, family, lam, head=None):
    pid = f"{family}_binary_lambda_{format(lam, '.8g').replace('.', 'p')}"
    directory = out / 'gate_training' / pid
    directory.mkdir(parents=True, exist_ok=True)
    completed = directory / 'complete.json'
    if completed.exists():
        require(read(completed)['signature'] == signature, 'Changed gate run')
        logging.info('Reuse completed %s', pid)
        return
    require(args.resume or not (directory / 'last.pt').exists(), 'Use --resume for interrupted gate training')
    seed_all(args.seed + 500)
    gate = BinaryGate(len(data.ids)).to(device)
    optimizer = torch.optim.Adam(gate.parameters(), lr=args.gate_lr, weight_decay=0.0)

    def forward(split, take):
        hard, _ = gate()
        if family == 'decoder':
            residual = head(data.raw.events(split, take, device), hard, fused=True)
        else:
            residual = profile_sum(csr_events(data.evidence[split], take, device), hard)
        return torch.tensor(data.explicit[split][take], dtype=torch.float32, device=device) + residual

    def evaluate():
        hard = (gate.logits.detach().sigmoid() > 0.5).cpu().numpy().astype(float)
        if family == 'profile':
            z = data.explicit['val'] + data.evidence['val'] @ hard
        else:
            z = np.empty(len(data.splits['val']), float)
            with torch.no_grad():
                for take in batches(len(z), args.batch_size):
                    z[take] = forward('val', take).cpu().numpy()
        m = metrics(data.splits['val'].labels, z)
        return (m, int(hard.sum()))
    history, start, stale, elapsed = ([], 0, 0, 0.0)
    if (directory / 'last.pt').exists():
        ck = load_pt(directory / 'last.pt')
        require(ck['signature'] == signature, 'Gate signature mismatch')
        gate.load_state_dict(ck['state'])
        optimizer.load_state_dict(ck['optimizer'])
        restore_rng(ck['rng'], device)
        history, start, stale, elapsed = (ck['history'], ck['epoch'], ck['stale'], ck['seconds'])
        best_state, best, best_epoch = (ck['best_state'], ck['best'], ck['best_epoch'])
    else:
        best, k = evaluate()
        best_state, best_epoch = (cpu_tree(gate.state_dict()), 0)
        history.append(dict(epoch=0, K=k, val_ap=best['ap'], val_auc=best['auc'], val_brier=best['brier']))
    if start < args.gate_epochs and stale < args.gate_patience:
        for epoch in range(start + 1, args.gate_epochs + 1):
            clock = time.monotonic()
            optimizer.zero_grad(set_to_none=True)
            total = 0.0
            split = data.splits['train']
            for step, take in enumerate(batches(len(split), args.batch_size)):
                z = forward('train', take)
                y = torch.tensor(split.labels[take], dtype=torch.float32, device=device)
                loss = F.binary_cross_entropy_with_logits(z, y, reduction='sum') / len(split)
                require(torch.isfinite(loss), 'Nonfinite gate BCE')
                loss.backward()
                total += float(loss.detach().cpu())
                if step and step % 200 == 0:
                    logging.info('%s epoch=%d microbatch=%d', pid, epoch, step)
            penalty = lam * gate.logits.sigmoid().mean()
            penalty.backward()
            require(gate.logits.grad is not None and torch.isfinite(gate.logits.grad).all(), 'Nonfinite gate gradient')
            optimizer.step()
            val, k = evaluate()
            if val['ap'] > best['ap'] + 1e-06:
                best, best_epoch, best_state, stale = (val, epoch, cpu_tree(gate.state_dict()), 0)
            else:
                stale += 1
            elapsed += time.monotonic() - clock
            history.append(dict(epoch=epoch, K=k, train_bce_before_update=total, sparsity_penalty_before_update=float(penalty.detach().cpu()), val_ap=val['ap'], val_auc=val['auc'], val_brier=val['brier'], seconds=elapsed))
            ck = dict(signature=signature, epoch=epoch, stale=stale, seconds=elapsed, state=cpu_tree(gate.state_dict()), optimizer=cpu_tree(optimizer.state_dict()), rng=rng_state(device), best_state=best_state, best=best, best_epoch=best_epoch, history=history)
            save_pt(directory / 'last.pt', ck)
            pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
            logging.info('%s epoch=%d/%d K=%d valAP=%.6f best=%d', pid, epoch, args.gate_epochs, k, val['ap'], best_epoch)
            release(device)
            if stale >= args.gate_patience:
                break
    gate.load_state_dict(best_state)
    hard = (gate.logits.detach().sigmoid() > 0.5).cpu().numpy().astype(float)
    logits = profile_logits(data, hard) if family == 'profile' else predict_head(data, head, args, device, torch.tensor(hard, dtype=torch.float32, device=device), fused=True)
    ref = f'{family}_full'
    record(out, data, pid, 'binary', family, hard, logits, reference_point=ref, gate_lambda=lam, best_epoch=best_epoch, gate_trainable_parameters=len(data.ids), gate_prob_mean=float(gate.logits.detach().sigmoid().mean().cpu()), seconds=elapsed)
    save_pt(directory / 'best.pt', dict(state=best_state, best_epoch=best_epoch, validation=best, feature_ids=torch.from_numpy(data.ids), signature=signature))
    write_json(completed, dict(signature=signature, validation_only=True))
    del gate, optimizer
    release(device)

def controls(args, data, out):
    points = data.proximal_points()
    budgets = sorted(set(args.budgets) | {int(p[2]['K']) for p in points})
    importance = np.asarray(abs(data.evidence['train']).mean(axis=0)).ravel()
    order = np.lexsort((data.ids, -importance))
    pd.DataFrame(dict(feature_id=data.ids[order], train_mean_abs_contribution=importance[order])).to_csv(out / 'train_contribution_ranking.tsv', sep='\t', index=False)
    fingerprint = {pid: provenance for pid, _, _, provenance in points}
    provenance_path = out / 'proximal_reference_manifest.json'
    if provenance_path.exists():
        require(read(provenance_path)['points'] == fingerprint, 'Proximal source changed; use a new E4 output root')
    write_json(provenance_path, dict(source=str(args.proximal_results), points=fingerprint, budgets=budgets))
    for pid, alpha, point, provenance in points:
        multipliers = np.divide(alpha, data.alpha, out=np.zeros_like(alpha), where=data.alpha > 0)
        require(not ((data.alpha == 0) & (alpha > 0)).any(), 'Nonzero amplitude for zero profile')
        name = 'profile_proximal_' + pid
        if not (out / f'points/{name}.json').exists():
            record(out, data, name, 'proximal', 'profile', multipliers, profile_logits(data, multipliers), reference_point='profile_full', converged=True, lambda_ratio=point['lambda_ratio'], provenance=provenance)
    for k in budgets:
        selected = np.zeros(len(data.ids))
        selected[order[:k]] = 1.0
        pid = f'profile_topk_K{k}'
        if not (out / f'points/{pid}.json').exists():
            record(out, data, pid, 'topk', 'profile', selected, profile_logits(data, selected), reference_point='profile_full', requested_K=k, selection_labels='TRAIN only')
        for repetition in range(args.random_repeats):
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, 804, repetition]))
            selected = np.zeros(len(data.ids))
            selected[rng.permutation(len(data.ids))[:k]] = 1.0
            pid = f'profile_random_r{repetition:02d}_K{k}'
            if not (out / f'points/{pid}.json').exists():
                record(out, data, pid, 'random', 'profile', selected, profile_logits(data, selected), reference_point='profile_full', requested_K=k, random_repeat=repetition, random_policy='uniform candidate permutation; repetitions are not performance-selected')

def refit_support(E, y, offset, ridge, max_iter, tol):
    if E.shape[1] == 0:
        return (np.zeros(0), dict(converged=True, iterations=0, kkt_max=0.0))

    def objective(x):
        z = offset + E @ x
        value = np.mean(np.logaddexp(0, z) - y * z) + 0.5 * ridge * (x @ x)
        gradient = np.asarray(E.T @ (expit(z) - y)).ravel() / len(y) + ridge * x
        require(np.isfinite(value) and np.isfinite(gradient).all(), 'Nonfinite refit objective')
        return (value, gradient)
    result = minimize(objective, np.ones(E.shape[1]), jac=True, bounds=[(0.0, None)] * E.shape[1], method='L-BFGS-B', options=dict(maxiter=max_iter, gtol=tol * 0.1, ftol=1e-14, maxls=50))
    value, gradient = objective(result.x)
    kkt = float(np.max(np.where(result.x > 1e-10, abs(gradient), np.maximum(-gradient, 0))))
    converged = kkt < tol
    require(np.isfinite(result.x).all(), 'Nonfinite refit coefficients')
    require(not 'ABNORMAL' in str(result.message).upper(), f'Refit line search failure: {result.message}')
    if not converged:
        logging.warning('Refit not converged: K=%d KKT=%.3g; excluded from frontier', E.shape[1], kkt)
    return (result.x, dict(converged=converged, iterations=int(result.nit), kkt_max=kkt, solver_success=bool(result.success), message=str(result.message), objective=value))

def refit(args, data, out):
    require((out / 'profile_complete.json').exists() and (out / 'controls_complete.json').exists(), 'Run profile and controls before refit')
    if not (out / 'points/profile_full_refit.json').exists():
        coef, diagnostic = refit_support(data.evidence['train'], data.splits['train'].labels, data.explicit['train'], args.refit_ridge, args.refit_max_iter, args.refit_tol)
        record(out, data, 'profile_full_refit', 'full', 'profile', coef, profile_logits(data, coef), reference_point='profile_full_refit', selected_support_K=len(data.ids), refit_ridge=args.refit_ridge, eligible=diagnostic['converged'], solver=diagnostic)
    paths = sorted((out / 'points').glob('*.json'))
    for path in paths:
        p = read(path)
        if p['family'] != 'profile' or p['method'] not in ['binary', 'proximal', 'topk', 'random']:
            continue
        pid = p['point_id'] + '_refit'
        if (out / f'points/{pid}.json').exists():
            continue
        original = load_pt(out / f"solutions/{p['point_id']}.pt")['coefficient'].numpy()
        support = np.flatnonzero(original > 1e-10)
        coef, diagnostics = refit_support(data.evidence['train'][:, support], data.splits['train'].labels, data.explicit['train'], args.refit_ridge, args.refit_max_iter, args.refit_tol)
        vector = np.zeros(len(data.ids))
        vector[support] = coef
        record(out, data, pid, p['method'] + '_refit', 'profile', vector, profile_logits(data, vector), reference_point='profile_full_refit', selected_support_K=len(support), source_point=p['point_id'], requested_K=p.get('requested_K'), random_repeat=p.get('random_repeat'), refit_ridge=args.refit_ridge, eligible=diagnostics['converged'], solver=diagnostics)

def main(argv=None):
    args = parse_args(argv)
    torch.set_num_threads(args.threads)
    device = device_for(args.device)
    out = args.output_root / f'{args.variant}_seed{args.seed}'
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / 'run.log')], force=True)
    with exclusive(out / 'process.lock'):
        logging.info('E4 phase=%s variant=%s seed=%d device=%s (Validation only)', args.phase, args.variant, args.seed, device)
        data = ExperimentData(args, with_raw=args.phase in ['check', 'decoder'])
        contract = dict(experiment='J_WARM_E4', training=training_contract(args), sources=data.fingerprint, code={p.name: digest(p) for p in [Path(__file__), HERE / 'e4_data.py', HERE / 'model_e4.py', HERE / 'experiment_helpers.py', HERE.parents[1] / 'train/train_jwarm.py', HERE.parents[1] / 'model/model_jwarm.py']}, explicit_source=data.source.name, gate='binary STE; logits=0 initialization', gate_penalty='lambda * mean(sigmoid(logits))', gate_optimizer_updates='one per full TRAIN epoch', candidate_count=len(data.ids), decoder_input='RAW latent, candidate-restricted; no P95', profile_input='J_WARM TRAIN P95 and post-joint frozen original-amplitude positional evidence', decoder_bias_included=True, explicit_frozen=True, gate_task_head_frozen=True, test_evaluated=False, stage_D=False, crop_mode='idr_coverage', context_margin=128)
        signature = identity(contract)
        config_path = out / 'config.json'
        if config_path.exists():
            require(read(config_path)['signature'] == signature, 'E4 inputs, code or hyperparameters differ from the saved run; use a new output root')
        write_json(config_path, dict(signature=signature, **contract))
        if args.phase == 'check':
            batch = data.raw.events('val', np.arange(2), torch.device('cpu'))
            write_json(out / 'input_check.json', dict(passed=True, rows={s: len(x) for s, x in data.splits.items()}, candidate_count=len(data.ids), raw_batch_nnz=len(batch[0]), decoder_vectors_shape=list(data.vectors.shape), no_Test_split_loaded=True, no_training_started=True))
            logging.info('E4 input check passed; no training started')
            return 0
        profile_reference(args, data, out)
        marker = out / (args.phase + '_complete.json')
        if marker.exists():
            require(read(marker)['signature'] == signature, 'Invalid completion marker')
            logging.info('Phase already complete; skipped')
            return 0
        if args.phase == 'decoder':
            decoder_stage(args, data, out, signature, device)
        elif args.phase == 'profile':
            for lam in args.gate_lambdas:
                fit_gate(args, data, out, signature, device, 'profile', lam)
        elif args.phase == 'controls':
            controls(args, data, out)
        elif args.phase == 'refit':
            refit(args, data, out)
        write_json(marker, dict(signature=signature, validation_only=True))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
