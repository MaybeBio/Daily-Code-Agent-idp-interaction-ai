import argparse, os, subprocess, sys
from pathlib import Path
p = argparse.ArgumentParser(description='Reproduce manuscript numerical panels from included result tables.')
p.add_argument('--figure', choices=['2', '3', '4', '5', '6', 'all'], default='all')
p.add_argument('--output', type=Path, default=Path('reproduce/output'))
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
e = os.environ.copy()
e['PRIORSLIM_PLOT_OUTPUT'] = str(a.output.resolve())
for f in '23456' if a.figure == 'all' else a.figure:
    subprocess.run([sys.executable, str(Path(__file__).resolve().parent / 'plots' / f'figure{f}.py')], env=e, check=True)
