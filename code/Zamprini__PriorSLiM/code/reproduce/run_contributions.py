from pathlib import Path
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser(description='Recompute latent contribution deletions on development folds.')
    p.add_argument('phase', choices=['check','discovery','followup','candidates'])
    p.add_argument('--data-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--fold', type=int, choices=range(1,11), default=6)
    args = p.parse_args()
    root, output = args.data_root.resolve(), args.output_dir.resolve()
    if output.is_relative_to(root): p.error('Output must be outside input data')
    os.environ['PRIORSLIM_DATA_ROOT'] = str(root)
    here = Path(__file__).resolve().parent
    sys.path[:0] = [str(here.parent/'src'),str(here/'ablation_helpers')]
    import _training
    sys.modules['train_jwarm'] = _training
    import latent_deletion as analysis
    analysis.OUT = output
    if args.phase == 'check':
        data = analysis.load_fold(args.fold)
        print(f"PASS: fold {args.fold}; {len(data['meta'])} validation windows; {len(data['ids'])} candidate features")
    elif args.phase == 'candidates':
        import latent_candidates as candidates
        candidates.OUT = output/'all_candidates'
        sys.argv = [sys.argv[0],'--phase','all']
        candidates.main()
    else:
        sys.argv = [sys.argv[0],'--phase',args.phase]
        analysis.main()

if __name__ == '__main__':
    main()
