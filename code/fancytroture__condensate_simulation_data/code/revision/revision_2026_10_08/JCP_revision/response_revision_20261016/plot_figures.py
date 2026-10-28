"""Response ratio figures: four observables, manuscript marker conventions."""
import os
os.environ['MPLCONFIGDIR']='/tmp/jcp_mpl_response_20261016'
from pathlib import Path
import csv,json,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];REV=ROOT/'revision_2026_10_08'
DATA=HERE.parent/'response_revision_20261015_v2'
OUT=HERE/'figures';OUT.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('mainstyle',ROOT/'project_figure/scripts/make_figure2_dynamics_story.py')
style=importlib.util.module_from_spec(spec);spec.loader.exec_module(style);style.setup_style()
SYSTEMS=['FUS','A1','TDP43'];NAMES=['FUS-LCD','A1-LCD','TDP-43-LCD']
METRICS=['rg','intra','inter','partner']
LABELS={'rg':r'$R_g$','intra':r'$N_{\mathrm{intra}}$','inter':r'$N_{\mathrm{inter}}$','partner':r'$S_{\mathrm{partner}}$'}
MARKERS={'rg':'o','intra':'D','inter':'^','partner':'s'}
FACES={'rg':'#4a4a4a','intra':'#4a4a4a','inter':'white','partner':'white'}
TAU=r'$\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$'
AREA=r'$I_W(T_{\mathrm{H}})/I_W(T_{\mathrm{L}})$'
VALUES=[]
LAYOUT=[]

def read(p):
    with p.open() as f:return list(csv.DictReader(f))

def legend(ax):
    handles=[Line2D([],[],marker=MARKERS[m],mfc=FACES[m],mec='#4a4a4a',mew=1.25,ls='none',color='#4a4a4a',ms=8.5,label=LABELS[m]) for m in METRICS]
    return ax.legend(handles=handles,loc='upper right',ncol=2,fontsize=15,
                     handletextpad=.4,columnspacing=1.1,borderaxespad=.6)

def axes_setup(ax,title,ylabel=TAU,letter=None,categorical=True,inside_title=False):
    ax.axhline(1,color='#777777',ls='--',lw=1,zorder=0)
    if inside_title:
        ax._inside_title=ax.text(.035,.965,title,transform=ax.transAxes,ha='left',va='top',
                                 fontsize=plt.rcParams['axes.titlesize'])
    else:
        ax.set_title(title,pad=13)
    ax.set_ylabel(ylabel)
    if categorical:ax.set_xticks(range(3),NAMES);ax.set_xlim(-.5,2.5)
    if letter:ax.text(-.13,1.04,letter,transform=ax.transAxes,fontsize=19,fontweight='bold')
    style.box_axis(ax)

def point(ax,x,row,metric,color='#4a4a4a',face=None,marker=None):
    val,lo,hi=(float(row[k]) for k in ['ratio','ci_low','ci_high'])
    assert np.isfinite([val,lo,hi]).all() and lo<=val<=hi,(row,metric)
    ax.errorbar(x,val,yerr=[[val-lo],[hi-val]],fmt=marker or MARKERS[metric],ms=8.5,
                mfc=face or FACES[metric],mec=color,mew=1.25,color=color,ecolor=color,
                elinewidth=1.05,capsize=5,capthick=1.05,zorder=3)

def save(fig,name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for i,ax in enumerate(fig.axes):
        legend=ax.get_legend()
        title=getattr(ax,'_inside_title',None)
        bounds={}
        for label,artist in [('legend',legend),('title',title)]:
            if artist is None:continue
            box=artist.get_window_extent(renderer).transformed(ax.transAxes.inverted())
            assert box.x0>=0 and box.y0>=0 and box.x1<=1 and box.y1<=1,(name,i,label,box)
            bounds[label]=[box.x0,box.y0,box.x1,box.y1]
        if legend is not None and title is not None:
            assert not legend.get_window_extent(renderer).overlaps(title.get_window_extent(renderer))
        if bounds:LAYOUT.append(dict(figure=name,panel=i,**bounds))
    for suffix in ['png','pdf']:fig.savefig(OUT/f'{name}.{suffix}',dpi=300,bbox_inches='tight')
    plt.close(fig)

def controls():
    data=read(REV/'paired_slab_temperature_jump/analysis_full20ns_2026_10_09/paired_dynamic_contrasts.csv')
    partners=read(HERE.parent/'response_revision_20261015/extended_partner/ratios.csv')
    fig,ax=plt.subplots(figsize=(9.2,4.8),layout='constrained')
    axes_setup(ax,'');ax.set_ylim(.25,1.35)
    for i,m in enumerate(METRICS):
        for x,s in enumerate(SYSTEMS):
            if m=='partner':row=dict(next(r for r in partners if r['system']==s and float(r['horizon_ns'])==17.5))
            else:
                k={'rg':'Rg','intra':'N_intra','inter':'N_inter'}[m]
                source=next(r for r in data if r['case']==s and r['window']=='full' and r['mode']=='original' and r['metric']==k and r['measure']=='threshold_time' and r['threshold']=='0.2')
                row={k:source[v] for k,v in [('ratio','high_over_low'),('ci_low','ratio_ci_low'),('ci_high','ratio_ci_high')]}
            point(ax,x+(i-1.5)*.14,row,m);VALUES.append({**row,'figure':'R1','metric':m,'system':s})
    legend(ax);save(fig,'Fig_R1_all_observable_controls')

def sensitivity():
    data=read(DATA/'sensitivity_ratios.csv')
    colors=['#2f7f73','#b54b5f','#436fa0'];markers=['o','s','^']
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.2),layout='constrained')
    for ax,m,letter in zip(axes.flat,METRICS,'ABCD'):
        axes_setup(ax,LABELS[m],letter=letter,categorical=False,inside_title=True)
        allhigh=[]
        for s,name,color,marker in zip(SYSTEMS,NAMES,colors,markers):
            rows=sorted([r for r in data if r['system']==s and r['metric']==m and r['kind']=='threshold'],key=lambda r:float(r['value']))
            xx=[float(r['value']) for r in rows];yy=[float(r['ratio']) for r in rows]
            ax.plot(xx,yy,color=color,lw=1.2,zorder=2)
            for x,row in zip(xx,rows):
                point(ax,x,row,m,color=color,face=(color if FACES[m]!='white' else 'white'))
                VALUES.append(dict(figure='R4',**row));allhigh.append(float(row['ci_high']))
        ax.set_xlabel(r'Retention level $\theta$' if m=='partner' else r'Correlation level $\theta$')
        ax.set_xticks([.4,.45,.5,.55,.6] if m=='partner' else [.1,.15,.2,.25,.3])
        ax.set_ylim(0,max(1.2,max(allhigh)*1.2))
        if letter=='A':ax.set_ylim(0,1.5)
    handles=[Line2D([],[],color=c,lw=2,label=n) for c,n in zip(colors,NAMES)]
    axes[0,0].legend(handles=handles,ncol=1,loc='upper right',fontsize=15,
                     bbox_to_anchor=(1,1.025),borderaxespad=.6)
    save(fig,'Fig_R4_all_observable_thresholds')
    fig,axes=plt.subplots(1,3,figsize=(15.,4.8),layout='constrained')
    for ax,s,name,letter in zip(axes,SYSTEMS,NAMES,'ABC'):
        axes_setup(ax,name,AREA,letter,categorical=False,inside_title=True)
        allhigh=[]
        for m in METRICS:
            rows=sorted([r for r in data if r['system']==s and r['metric']==m and r['kind']=='integral'],key=lambda r:float(r['value']))
            xs=[float(r['value']) for r in rows];ys=[float(r['ratio']) for r in rows]
            ax.plot(xs,ys,lw=1,color='#777777')
            for x,row in zip(xs,rows):
                point(ax,x,row,m);VALUES.append(dict(figure='R5',**row));allhigh.append(float(row['ci_high']))
        ax.set_xticks([10,20,30]);ax.set_xlim(7,33);ax.set_xlabel(r'Integration window $W$ (ns)')
        ax.set_ylim(0,max(1.3,max(allhigh)*1.25))
        if letter=='A':ax.set_ylim(0,1.95)
    legend(axes[0]);save(fig,'Fig_R5_all_observable_integrals')

def residence():
    data=read(DATA/'residence/new_ratios.csv')+read(REV/'R2C5_chain_selection/central60_2026_10_10/temperature_ratios.csv')
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.8),layout='constrained')
    for ax,group,title,letter in zip(axes,['original_selection','short_matched_original'],['Original selection','100 ns residence'],'AB'):
        axes_setup(ax,title,letter=letter,inside_title=True);allhigh=[]
        for i,m in enumerate(METRICS):
            for x,s in enumerate(SYSTEMS):
                row=next(r for r in data if r['system']==s and r['metric']==m and r['group']==group and r['quantity']=='tau')
                point(ax,x+(i-1.5)*.14,row,m)
                VALUES.append(dict(figure='R7',**row));allhigh.append(float(row['ci_high']))
        ax.set_ylim(0,max(allhigh)*(1.4 if letter=='B' else 1.17))
    legend(axes[1]);save(fig,'Fig_R7_all_observable_residence')

def mobility():
    data=read(REV/'R2C5_chain_selection/central60_2026_10_10/mobility_summary.csv')
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.8),layout='constrained')
    temperatures=[('FUS','280','325'),('A1','230','260'),('TDP43','260','300')]
    for j,selected in enumerate(['original_window_matched','original_ids_bulk_origins']):
        ax=axes[j]
        ax.text(-.13,1.04,'AB'[j],transform=ax.transAxes,fontsize=19,fontweight='bold')
        style.box_axis(ax)
        for role,offset in [('Tl',-.12),('Th',.12)]:
            for i,(system,low,high) in enumerate(temperatures):
                temp=low if role=='Tl' else high
                row=next(r for r in data if r['case']==f'{system}_{temp}' and r['lag_ns']=='10'
                         and r['coordinates']=='xy' and r['selected']==selected)
                val,lo,hi=(float(row[k]) for k in ['ratio','ci_low','ci_high'])
                ax.errorbar(i+offset,val,yerr=[[max(0,val-lo)],[max(0,hi-val)]],fmt='o',
                            color=style.TEMP_COLORS[role],capsize=4,ms=7,lw=1.5)
        ax.axhline(1,color='.45',ls=':',lw=1.1)
        ax.set_xticks(range(3),NAMES);ax.set_xlim(-.5,2.5)
        ax.set_ylabel(r'$M_{xy}^{\mathrm{selected}}(10)/M_{xy}^{\mathrm{all}}(10)$')
        ax.set_ylim(.55,1.3);ax.yaxis.set_major_locator(MultipleLocator(.2))
    axes[0].legend(handles=[Line2D([],[],color=style.TEMP_COLORS[r],marker='o',ls='none',
                                  label=style.temp_label(r)) for r in ['Tl','Th']],
                   loc='upper left',ncol=2)
    save(fig,'Fig_R6_mobility')

if __name__=='__main__':
    controls();sensitivity();residence();mobility()
    (HERE/'plotted_values.json').write_text(json.dumps(VALUES,indent=2)+'\n')
    (HERE/'layout_bounds.json').write_text(json.dumps(LAYOUT,indent=2)+'\n')
    for fig in ['R1','R4','R5','R7']:
        assert {r['metric'] for r in VALUES if r['figure']==fig}==set(METRICS)
    print('All four dynamics figures contain all four observables and 95% intervals.')
