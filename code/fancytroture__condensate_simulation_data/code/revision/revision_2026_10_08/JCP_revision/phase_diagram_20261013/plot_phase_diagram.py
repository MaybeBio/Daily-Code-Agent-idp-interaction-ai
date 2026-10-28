#!/usr/bin/env python3
"""Plot existing coexistence densities and their existing density-gap fit."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-jcp-phase')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / 'critical_temperature_slab_simulations/density_3to7us_20261012/si_figure'

def read(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))

densities = read('critical_temperature_input_densities.csv')
fits = {r['system']: r for r in read('critical_temperature_fit.csv')}
STYLE = HERE.parents[2] / 'project_figure/scripts/make_figure2_dynamics_story.py'
spec = importlib.util.spec_from_file_location('maintext_style', STYLE)
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
style.setup_style()
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
fig.subplots_adjust(left=.075, right=.985, bottom=.20, top=.77, wspace=.30)
blue = '#333333'
verification = {}
for j, (system, title) in enumerate([('FUS','FUS-LCD'),('A1','A1-LCD'),('TDP43','TDP-43-LCD')]):
    rows = sorted((r for r in densities if r['system']==system),
                  key=lambda r:float(r['temperature_K']))
    temp = np.array([float(r['temperature_K']) for r in rows])
    dense = np.array([float(r['dense_g_ml']) for r in rows])
    dilute = np.array([float(r['dilute_g_ml']) for r in rows])
    tc, amp, beta = (float(fits[system][k]) for k in ('Tc_K','amplitude','beta'))
    gap = dense - dilute
    rmse = float(np.sqrt(np.mean((gap-amp*(tc-temp)**beta)**2)))
    assert abs(rmse-float(fits[system]['rmse_g_ml'])) < 1e-10
    verification[system] = dict(temperatures_K=temp.tolist(), Tc_K=tc,
                                rmse_g_ml=rmse, density_pairs=len(rows))
    ax = axes[j]
    diameter = (dense+dilute)/2
    slope, rho_c = np.polyfit(tc-temp,diameter,1)
    verification[system]['critical_density_g_ml'] = float(rho_c)
    verification[system]['diameter_slope_g_ml_K'] = float(slope)
    # Solid curves interpolate the measured density branches. No new trajectory
    # data or dense-minus-dilute critical fit is introduced.
    t_obs = np.linspace(temp[0],temp[-1],250)
    for branch in (dilute,dense):
        ax.plot(PchipInterpolator(temp,branch)(t_obs),t_obs,color=blue,lw=2.3)
    point_colors = [style.TEMP_COLORS['Tl'],style.TEMP_COLORS['Th'],blue]
    for index, color in enumerate(point_colors):
        ax.plot(dilute[index],temp[index],'o',color=color,mfc='white',ms=7.5,zorder=3)
        ax.plot(dense[index],temp[index],'o',color=color,mfc=color,ms=7.5,zorder=3)
    for index, label, offset, align in [(0,r'$T_{\mathrm{L}}$',(-12,8),'right'),
                                        (1,r'$T_{\mathrm{H}}$',(10,5),'left')]:
        ax.annotate(label,(dense[index],temp[index]),xytext=offset,
                    textcoords='offset points',ha=align,va='bottom',
                    color=point_colors[index],fontsize=14.8)
    # The dashed closure is a scaling guide anchored to the highest-T pair,
    # using the existing Tc/beta and a rectilinear-diameter estimate of rho_c.
    u = np.linspace(0,1,400)
    reduced_t = u**(1/beta)
    t_cap = tc-(tc-temp[-1])*reduced_t
    center_cap = rho_c+(diameter[-1]-rho_c)*reduced_t
    gap_cap = gap[-1]*u
    for branch in (center_cap-.5*gap_cap,center_cap+.5*gap_cap):
        assert branch.min() >= -1e-12
        ax.plot(branch,t_cap,color=blue,lw=2.3,ls='--')
    ax.plot(rho_c,tc,'D',color=blue,ms=6.5,zorder=4)
    ax.text(rho_c+.018,tc+3.4,rf'$T_{{\mathrm{{c}}}}\approx {tc:.0f}\,\mathrm{{K}}$',
            ha='left',va='center',fontsize=14.8)
    ax.set(xlim=(-.028,.66),ylim=(temp[0]-6,tc+9),
           xlabel=r'$\rho$ (g/mL)',ylabel=r'$T$ (K)')
    ax.set_xticks([0,.2,.4,.6])
    ax.set_title(title,pad=12)
    ax.text(-.13,1.055,'ABC'[j],transform=ax.transAxes,weight='bold',fontsize=20)
    style.box_axis(ax)

fig.legend(handles=[Line2D([],[],marker='o',ls='none',color=blue,mfc='white',label='Dilute phase'),
                    Line2D([],[],marker='o',ls='none',color=blue,label='Dense phase')],
           loc='upper center',bbox_to_anchor=(.5,1.0),ncol=2,frameon=False,
           handletextpad=.5,columnspacing=2)
for suffix in ('pdf','png'):
    fig.savefig(HERE / ('Fig_SI_phase_diagram.'+suffix),facecolor='white',dpi=300)
plt.close(fig)
verification['input_sha256'] = {name:hashlib.sha256((DATA/name).read_bytes()).hexdigest()
    for name in ('critical_temperature_input_densities.csv','critical_temperature_fit.csv')}
verification['style_source'] = str(STYLE)
verification['style'] = {key:matplotlib.rcParams[key] for key in
    ('font.family','axes.labelsize','axes.titlesize','xtick.labelsize','ytick.labelsize',
     'legend.fontsize','axes.linewidth','pdf.fonttype')}
(HERE/'phase_data_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
print(json.dumps(verification,indent=2))
