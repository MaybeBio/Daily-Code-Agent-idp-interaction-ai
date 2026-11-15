from pathlib import Path
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE = Path(__file__).resolve().parents[1]
DATA = BASE / 'data'
OUT = Path(os.environ.get('PRIORSLIM_PLOT_OUTPUT', str(BASE / 'output')))
OUT.mkdir(parents=True, exist_ok=True)

def save(fig, figure, panel):
    for ext in ('svg', 'png'):
        fig.savefig(OUT / f'fig{figure}_{panel}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
