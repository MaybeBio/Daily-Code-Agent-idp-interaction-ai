from pathlib import Path
import csv,json,hashlib,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent.parent; W=Path(__file__).resolve().parent; R=P.parent
OLD=R/'paired_slab_temperature_jump/analysis_full20ns_2026_10_09'
CASES=[('FUS','FUS-LCD',280,325),('A1','A1-LCD',230,260),('TDP43','TDP-43-LCD',260,300)]
weights=np.load(OLD/'bootstrap_pair_weights.npz')['weights']
with (OLD/'per_snapshot_time_evolution.csv').open() as h: old=list(csv.DictReader(h))
cache=W/'density_summary.npz'
if cache.exists():
 data=dict(np.load(cache))
else:
 data={}; hashes={}
 for case,title,low,high in CASES:
  for temp in (low,high):
   summaries=[]
   for snap in range(1,9):
    root=R/f'{case}_slab_temperature_jump/runs/snapshot{snap:02d}_{temp}K'
    branches=[];times=[]
    for part in sorted(root.glob('part[0-9][0-9][0-9][0-9]')):
     path=part/'density.npz';sha=hashlib.sha256(path.read_bytes()).hexdigest()
     assert sha==json.loads((part/'complete.json').read_text())['files']['density.npz']
     hashes[str(path)]=sha
     with np.load(path) as d:
      z=(d['edges_nm'][:-1]+d['edges_nm'][1:])/2;prof=d['profiles']
      dense=prof[:,abs(z)<5].mean(1);dilute=prof[:,abs(z)>80].mean(1)
      np.testing.assert_allclose(dense,d['records'][:,2],atol=1e-12)
      branches.append(np.stack([dense,dilute],axis=1));times.append(d['records'][:,0])
    branch=np.concatenate(branches);time=np.concatenate(times)
    np.testing.assert_allclose(time,np.arange(1,2001)/100,atol=1e-5)
    time=np.arange(1,2001)/100  # Match the existing exact saved-frame binning.
    summary=np.stack([branch[(time>max(.1,start)+1e-8)&(time<=start+.5+1e-8)].mean(0) for start in np.arange(0,20,.5)])
    existing=sorted((r for r in old if r['metric']=='density' and r['case']==case and int(r['temperature_K'])==temp and int(r['snapshot'])==snap),key=lambda r:float(r['time_ns']))
    np.testing.assert_allclose(summary[:,0],[float(r['value']) for r in existing],rtol=0,atol=1e-12)
    summaries.append(summary)
   data[f'{case}_{temp}']=np.stack(summaries)
  print('Read',case,flush=True)
 np.savez_compressed(cache,**data)
 (W/'profile_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
spec=importlib.util.spec_from_file_location('style',P.parents[1]/'project_figure/scripts/make_figure2_dynamics_story.py');style=importlib.util.module_from_spec(spec);spec.loader.exec_module(style);style.setup_style()
fig,axes=plt.subplots(2,3,figsize=(15.5,7.0),sharex=True)
times=np.arange(.25,20,.5);report={}
for col,(case,title,low,high) in enumerate(CASES):
 for row,phase in enumerate(('dense','dilute')):
  ax=axes[row,col];ymax=0
  for temp,role in ((low,'Tl'),(high,'Th')):
   values=data[f'{case}_{temp}'][:,:,row];mean=values.mean(0);ci=np.quantile(weights@values,[.025,.975],axis=0)
   color=style.TEMP_COLORS[role];ax.plot(times,mean,color=color,lw=2,label=style.temp_label(role));ax.fill_between(times,*ci,color=color,alpha=.16,linewidth=0)
   ymax=max(ymax,float(ci[1].max()))
   report[f'{case}_{temp}_{phase}']={'min':float(mean.min()),'max':float(mean.max()),'mean':float(mean.mean()),'bootstrap_ci_max':float(ci[1].max())}
  ax.set(xlim=(0,20),ylabel=rf'$\rho_{{\mathrm{{{phase}}}}}$ (g/mL)');ax.set_xticks([0,5,10,15,20]);style.box_axis(ax)
  ax.text(.025,.965,chr(65+row*3+col),transform=ax.transAxes,fontsize=20,fontweight='bold',va='top')
  if row==0:ax.set_title(title,pad=12)
  else:
   ax.set_xlabel('Time after branching (ns)')
   ax.set_ylim(0,1.15*ymax if ymax else 1e-5)
   ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
  if row==0 and col==0:ax.legend(loc='upper right',ncol=2,handlelength=1.4,columnspacing=.7)
fig.subplots_adjust(left=.085,right=.985,bottom=.13,top=.91,wspace=.40,hspace=.30)
for ext in ('pdf','png'):fig.savefig(W/'figures'/f'Fig_control_density_only.{ext}',dpi=220)
plt.close(fig)
(W/'density_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
