from pathlib import Path
import argparse
import importlib
import os
import sys


def main():
    parser = argparse.ArgumentParser(description='Run fold-6 training and feature-selection experiments.')
    parser.add_argument('phase', choices=['check', 'prepare', 'smoke', 'train', 'compress', 'refine', 'decoder', 'profile', 'controls', 'refit', 'analyze', 'analysis-check'])
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--arm', choices=['A','B','C','D','B1','B2','B3','CE','S','J_SCRATCH','H3','AROMATIC_ADD','AROMATIC_SPLIT'], default='C')
    parser.add_argument('--seed', type=int, choices=range(2026,2031), default=2026)
    parser.add_argument('--device', choices=['cpu','mps','cuda','auto'], default='cpu')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--decoder-checkpoint', type=Path, help='Official MotifAE k3 weights for decoder-based controls')
    args = parser.parse_args()
    root = args.data_root.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    if output.is_relative_to(root):
        parser.error('--output-dir must be outside --data-root')
    os.environ['PRIORSLIM_DATA_ROOT'] = str(root)
    base = Path(__file__).resolve().parent
    sys.path[:0] = [str(base.parent/'src'), str(base/'ablation_helpers')]
    import _training, model, train
    sys.modules.update(train_jwarm=_training, model_jwarm=model, train_cv10=train)
    if args.phase in ['analyze','analysis-check']:
        import analyze_fold06 as analysis
        analysis.RESULTS = output/'training'
        analysis.OUT = output/'analysis'
        analysis.dep.E1 = analysis.RESULTS
        analysis.dep.self_test()
        if args.phase == 'analysis-check':
            print('PASS: group deletion numerical self-test')
        else:
            analysis.main()
        return
    if args.phase in ['decoder','profile','controls','refit']:
        import run_e4_fold06 as wrapper
        engine = wrapper.engine
        engine.ExperimentData = wrapper.FoldE4Data
        command = [args.phase,'--seed',str(args.seed),'--variant','C','--device',args.device,
                   '--project-root',str(root),'--e1-results',str(output/'training'),
                   '--proximal-results',str(output/'training'),'--output-root',str(output/'selection'),
                   '--budgets','64','128','256','512','1024','2048','4096','5829']
        if args.decoder_checkpoint:
            command.extend(['--decoder-checkpoint', str(args.decoder_checkpoint.expanduser().resolve())])
        if args.resume: command.append('--resume')
        engine.main(command)
        return
    import train_fold06 as experiment
    command = ['prepare' if args.phase=='check' else args.phase,'--arm',args.arm,'--seed',str(args.seed),
               '--device',args.device,'--output-root',str(output/'training')]
    if args.resume: command.append('--resume')
    if args.phase=='check':
        options = experiment.arguments(command)
        train.verify_fold_contract(options)
        for relative in ['train/combined_tensor.npz','train/combined_feature_manifest.json',
                         'feature_cache/motifae_sparse_windows_idr_coverage/motifae_sparse_config.json',
                         'train/outputs_jwarm_cv10/seed42/fold_06/_prepared/fold_training_assets.pt']:
            if not (root/relative).is_file(): raise FileNotFoundError(root/relative)
        print('PASS: fold-6 split hashes and required asset paths')
    else:
        experiment.main(command)

if __name__ == '__main__':
    main()
