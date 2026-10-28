"""Retain the requested selected-versus-all mobility comparison as Figure R6."""
import os
os.environ['MPLCONFIGDIR']='/tmp/jcp_mpl_remove_R6B_20261016'
from pathlib import Path
import csv
import json
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
REV=ROOT/'revision_2026_10_08'
OUT=HERE/'figures'
OUT.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('mainstyle',ROOT/'project_figure/scripts/make_figure2_dynamics_story.py')
style=importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
style.setup_style()
with (REV/'R2C5_chain_selection/central60_2026_10_10/mobility_summary.csv').open() as f:
    data=list(csv.DictReader(f))
temperatures=[('FUS','280','325'),('A1','230','260'),('TDP43','260','300')]
fig,ax=plt.subplots(figsize=(6.4,4.8),layout='constrained')
style.box_axis(ax)
rows=[]
for role,offset in [('Tl',-.12),('Th',.12)]:
    for i,(system,low,high) in enumerate(temperatures):
        temp=low if role=='Tl' else high
        row=next(r for r in data if r['case']==f'{system}_{temp}' and r['lag_ns']=='10'
                 and r['coordinates']=='xy' and r['selected']=='original_window_matched')
        val,lo,hi=(float(row[k]) for k in ['ratio','ci_low','ci_high'])
        assert lo<=val<=hi
        rows.append(row)
        ax.errorbar(i+offset,val,yerr=[[val-lo],[hi-val]],fmt='o',
                    color=style.TEMP_COLORS[role],capsize=4,ms=7,lw=1.5)
ax.axhline(1,color='.45',ls=':',lw=1.1)
ax.set_xticks(range(3),['FUS-LCD','A1-LCD','TDP-43-LCD'])
ax.set_xlim(-.5,2.5)
ax.set_ylabel(r'$M_{xy}^{\mathrm{selected}}(10)/M_{xy}^{\mathrm{all}}(10)$')
ax.set_ylim(.55,1.3)
ax.yaxis.set_major_locator(MultipleLocator(.2))
legend=ax.legend(handles=[Line2D([],[],color=style.TEMP_COLORS[r],marker='o',ls='none',
                               label=style.temp_label(r)) for r in ['Tl','Th']],
                 loc='upper left',ncol=2)
fig.canvas.draw()
box=legend.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transAxes.inverted())
assert box.x0>=0 and box.y0>=0 and box.x1<=1 and box.y1<=1
for suffix in ['png','pdf']:
    fig.savefig(OUT/f'Fig_R6_mobility.{suffix}',dpi=300,bbox_inches='tight')
plt.close(fig)
(HERE/'retained_mobility_values.json').write_text(json.dumps(rows,indent=2)+'\n')
print('Saved single-panel R6 with the original six values and confidence intervals.')
