#!/usr/bin/env python3
"""Compare original, index-only, and fully realigned TDP-43 state distances."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sys.path.insert(0,str(OUT.parent/'followup'))
import analyze_followup as previous
import rebuild_contacts as rebuild
s=previous.s
SAM=previous.SAM
NBOOT=4000
SEED=20261009


def write(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader();w.writerows(rows)


def permute_statistics(c,m,perm):
    return (np.stack([x[...,p,:][...,:,p] for x,p in zip(c,perm)]),
            np.stack([x[...,p] for x,p in zip(m,perm)]))


def sample_condition(family,role,rng):
    tag=f'TDP43_{family}_{s.TEMPS["TDP43"][role]}'
    with np.load(OUT/tag/'fixed_alignment.npz') as z:
        perm=z['new_to_old'];fixed_mu=z['avg_mu'];attempts=z['attempts']
    with np.load(SAM/f'{tag}_sufficient_statistics.npz') as z:
        assert np.array_equal(attempts,z['attempts'])
        c=z['counts'].sum(2);m=z['masses'].sum(2)
        assert int(z['lag'])==200 and float(z['dt'])==.01
    full_c,full_m=c.sum(1),m.sum(1)
    fc,fm=permute_statistics(full_c,full_m,perm)
    np.testing.assert_allclose(previous.pm(fm).mean(0),fixed_mu,atol=1e-6)
    original=previous.old.load_full(previous.FULL,'TDP43',family,role)
    np.testing.assert_allclose(previous.transition(full_c).mean(0),original['t'],atol=1e-12)
    result={v:{} for v in ['original','fixed']}
    nm=len(perm)
    weights=rng.multinomial(nm,np.full(nm,1/nm),NBOOT)/nm
    for version,ct,mt in [('original',full_c,full_m),('fixed',fc,fm)]:
        t,p=previous.transition(ct),previous.pm(mt)
        result[version]['model_seeds']=(np.concatenate([t.mean(0)[None],np.einsum('bm,mij->bij',weights,t)]),
                                       np.concatenate([p.mean(0)[None],weights@p]))
    blocks=[0] if family=='free' else [10,50,100]
    for block in blocks:
        if family=='condensate':
            with np.load(SAM/f'{tag}_absolute_statistics.npz') as z:
                assert len(z['excluded_trajectories'])==0
                factor=block//10
                c=z['counts'].reshape(nm,-1,factor,5,5).sum(2)
                m=z['masses'].reshape(nm,-1,factor,5).sum(2)
        nu=c.shape[1]
        w=rng.multinomial(nu,np.full(nu,1/nu),NBOOT)
        key='sampling' if block in [0,50] else f'sampling_{block}ns'
        corrected=permute_statistics(c,m,perm)
        for version,(ct,mt) in [('original',(c,m)),('fixed',corrected)]:
            bt=previous.transition(np.einsum('bu,muij->bmij',w,ct,optimize=True)).mean(1)
            bp=previous.pm(np.einsum('bu,mui->bmi',w,mt,optimize=True)).mean(1)
            result[version][key]=(np.concatenate([previous.transition(ct.sum(1)).mean(0)[None],bt]),
                                  np.concatenate([previous.pm(mt.sum(1)).mean(0)[None],bp]))
    return result


def contact_pair(family,variant):
    pairs=[]
    for role in s.ROLES:
        tag=f'TDP43_{family}_{s.TEMPS["TDP43"][role]}'
        if variant=='fully_fixed':
            with np.load(OUT/tag/'fixed_alignment.npz') as z:
                pairs.append((z['avg_contacts'],z['avg_mu']))
        else:
            d=s.model_dir('TDP43',family,role)
            a=np.load(d/'avg_contacts.npy').astype(float)
            if variant=='index_only' and role=='Tl':a=previous.indexing_diagnostic(a)
            pairs.append((a,np.load(d/'avg_mu.npy').astype(float)))
    return pairs


def match(family,variant,cutoff=.90):
    pairs=contact_pair(family,variant)
    n=min(pairs[0][0].shape[-1],pairs[1][0].shape[-1])
    f=[previous.contact_features(a[:,:n,:n]) for a,p in pairs]
    corr=f[1]@f[0].T
    kept=[np.where(p>=.10)[0] for a,p in pairs]
    rr,cc=linear_sum_assignment(1-corr[np.ix_(kept[1],kept[0])])
    accepted=[(int(kept[0][j]),int(kept[1][i])) for i,j in zip(rr,cc)
              if corr[kept[1][i],kept[0][j]]>=cutoff]
    # Display identities anchored to low-temperature local state rank.
    mapping=[{l:l for l in kept[0]},{}]
    assigned={h:l for l,h in accepted}
    next_id=5
    for h in kept[1]:
        if int(h) in assigned:mapping[1][assigned[int(h)]]=int(h)
        else:mapping[1][next_id]=int(h);next_id+=1
    return mapping,corr,accepted,pairs


def all_metrics(arrays,mapping):
    out,common,ids,ts,ps=previous.metrics(arrays,mapping)
    w=.5*(previous.pm(ps[0])+previous.pm(ps[1]))
    exit_c=[1-t.sum(-1) for t in ts]
    extra=np.sum(w*(exit_c[1]-exit_c[0])**2,-1)
    out['D_common_plus_exit']=np.sqrt(out['D_common_pooled_weights']**2+extra)
    cond=[t/t.sum(-1,keepdims=True) for t in ts]
    out['D_common_conditional']=np.sqrt(np.sum(w*np.sum((cond[1]-cond[0])**2,-1),-1))
    out['common_population_TL']=ps[0].sum(-1)
    out['common_population_TH']=ps[1].sum(-1)
    return out,common,ts,ps,exit_c


def analyze(data):
    rows,contrasts,matching,entries,populations,matrices=[],[],[],[],[],[]
    draws={}
    for variant in ['original','index_only','fully_fixed']:
        for cutoff in ([.90] if variant!='fully_fixed' else [.80,.85,.90,.95]):
            for family in s.FAMILIES:
                mapping,corr,accepted,pairs=match(family,variant,cutoff)
                if variant=='original':
                    reference=previous.maps('TDP43',family)
                    a={(reference[0][x],reference[1][x]) for x in set(reference[0])&set(reference[1])}
                    assert set(accepted)==a,('Cannot reproduce original cross-temperature mapping',family,accepted,a)
                    mapping=reference
                for hi in range(5):
                    for lo in range(5):
                        matching.append(dict(variant=variant,family=family,cutoff=cutoff,local_TL=lo,local_TH=hi,
                            contact_correlation=float(corr[hi,lo]),accepted=(lo,hi) in accepted,
                            selected_TL=bool(pairs[0][1][lo]>=.10),selected_TH=bool(pairs[1][1][hi]>=.10)))
                version='fixed' if variant=='fully_fixed' else 'original'
                for source in data[family,'Tl'][version]:
                    arrays=[data[family,r][version][source] for r in s.ROLES]
                    out,common,ts,ps,exit_c=all_metrics(arrays,mapping)
                    for metric,vals in out.items():
                        key=(variant,cutoff,family,source,metric)
                        draws[key]=vals
                        rows.append(dict(variant=variant,cutoff=cutoff,family=family,uncertainty=source,
                            metric=metric,common_pairs=str(accepted),n_common=len(common),**previous.summary(vals)))
                    if variant=='fully_fixed' and cutoff==.90:
                        for i,state in enumerate(common):
                            base=dict(family=family,uncertainty=source,source_state=state+1,lag_ns=2.)
                            populations.append(dict(**base,population_TL=float(ps[0][0,i]),population_TH=float(ps[1][0,i]),
                                **previous.summary(ps[1][:,i]-ps[0][:,i])))
                            for j,target in enumerate(common):
                                entries.append(dict(**base,target_state=target+1,T_TL=float(ts[0][0,i,j]),T_TH=float(ts[1][0,i,j]),
                                    **previous.summary(ts[1][:,i,j]-ts[0][:,i,j])))
                        if source=='model_seeds':
                            for role,(t,p) in zip(s.ROLES,arrays):
                                for i in range(5):
                                    for j in range(5):
                                        matrices.append(dict(family=family,role=role,local_source=i,local_target=j,
                                                             population=float(p[0,i]),transition_probability=float(t[0,i,j])))
            for source in ['sampling','model_seeds','sampling_10ns','sampling_100ns']:
                free_source=source if source in ['sampling','model_seeds'] else 'sampling'
                for metric in out:
                    diff=draws[variant,cutoff,'condensate',source,metric]-draws[variant,cutoff,'free',free_source,metric]
                    contrasts.append(dict(variant=variant,cutoff=cutoff,uncertainty=source,metric=metric,**previous.summary(diff)))
    write('distance_results.csv',rows);write('condensate_minus_free.csv',contrasts)
    write('contact_matching.csv',matching);write('fixed_population_changes.csv',populations)
    write('fixed_transition_entry_changes.csv',entries);write('fixed_full_local_matrices.csv',matrices)
    np.savez_compressed(OUT/'metric_bootstrap_draws.npz',**{'__'.join(map(str,k)):v for k,v in draws.items()})
    return rows,contrasts


def plot_comparison(rows,contrasts):
    previous.style()
    fig,axs=plt.subplots(1,3,figsize=(14,5.4))
    variants=['original','index_only','fully_fixed']
    labels=['Original\nmapping','Index only','Index +\nseed alignment']
    for ax,metric,title in zip(axs[:2],['D_common_pooled_weights','D_common_uniform_weights'],
                              ['Population-weighted distance','Uniform-weight distance']):
        for x,variant in enumerate(variants):
            for shift,family in [(-.16,'free'),(.16,'condensate')]:
                for offset,source,color in [(-.045,'sampling','#303030'),(.045,'model_seeds','#397eb1')]:
                    r=previous.old.one(rows,variant=variant,cutoff=.90,family=family,uncertainty=source,metric=metric)
                    xx=x+shift+offset
                    ax.plot([xx,xx],[r['ci_low'],r['ci_high']],color=color,lw=1.5)
                    ax.plot(xx,r['estimate'],'o',color=color,mfc='white' if family=='free' else color,ms=6)
        ax.set_title(title,pad=12);ax.set_ylim(0,.25);ax.set_ylabel('Common-state distance')
    ax=axs[2]
    for x,variant in enumerate(variants):
        for offset,source,color in [(-.07,'sampling','#303030'),(.07,'model_seeds','#397eb1')]:
            r=previous.old.one(contrasts,variant=variant,cutoff=.90,uncertainty=source,metric='D_common_pooled_weights')
            ax.plot([x+offset]*2,[r['ci_low'],r['ci_high']],color=color,lw=1.5)
            ax.plot(x+offset,r['estimate'],'o',color=color,ms=6)
    ax.axhline(0,color='.5',ls=':',lw=1);ax.set_ylim(-.10,.05)
    ax.set_ylabel(r'$D_{\mathrm{condensate}}-D_{\mathrm{free}}$');ax.set_title('Weighted-distance contrast',pad=12)
    for ax in axs:ax.set_xticks(range(3),labels);ax.set_xlim(-.5,2.5)
    handles=[Line2D([],[],color='black',marker='o',mfc='white',ls='none',label='Free chain'),
             Line2D([],[],color='black',marker='o',ls='none',label='Condensate'),
             Line2D([],[],color='#303030',lw=2,label='Sampling'),Line2D([],[],color='#397eb1',lw=2,label='Model seeds')]
    fig.suptitle('TDP-43: contact-index correction',fontsize=20,y=.99)
    fig.legend(handles=handles,loc='upper center',ncol=4,bbox_to_anchor=(.5,.92),fontsize=12)
    fig.subplots_adjust(left=.065,right=.98,top=.75,bottom=.20,wspace=.32)
    fig.text(.5,.035,'All distances use unfiltered transitions. Matching: population >= 0.10; correlation >= 0.90. Separate 95% intervals.',
             ha='center',fontsize=10)
    for ext in ['png','pdf']:fig.savefig(OUT/f'TDP43_before_after.{ext}',dpi=200)
    plt.close(fig)


def plot_matching():
    previous.style()
    fig,axes=plt.subplots(2,2,figsize=(9.4,8.4))
    for row,variant in enumerate(['original','fully_fixed']):
        for col,family in enumerate(s.FAMILIES):
            ax=axes[row,col]
            mapping,corr,accepted,pairs=match(family,variant,.90)
            ids=[np.where(p>=.10)[0] for a,p in pairs]
            values=corr[np.ix_(ids[1],ids[0])]
            im=ax.imshow(values,vmin=.6,vmax=1.,cmap='YlGnBu',aspect='equal')
            for i,hi in enumerate(ids[1]):
                for j,lo in enumerate(ids[0]):
                    chosen=(int(lo),int(hi)) in accepted
                    ax.text(j,i,f'{values[i,j]:.3f}',ha='center',va='center',fontsize=13,
                            color='white' if values[i,j]>.85 else 'black',fontweight='bold' if chosen else 'normal',
                            bbox=dict(boxstyle='square,pad=.18',fc='none',ec='white',lw=1.4) if chosen else None)
            ax.set_xticks(range(len(ids[0])),[f'{i+1}' for i in ids[0]])
            ax.set_yticks(range(len(ids[1])),[f'{i+1}' for i in ids[1]])
            ax.set_xlabel('Low-temperature local state')
            ax.set_ylabel(('Original' if row==0 else 'Corrected')+'\nHigh-temperature local state')
            if row==0:ax.set_title('Free chain' if family=='free' else 'Condensate',pad=12)
    fig.subplots_adjust(left=.11,right=.86,bottom=.10,top=.89,wspace=.38,hspace=.32)
    cb=fig.add_axes([.91,.25,.023,.5]);fig.colorbar(im,cax=cb)
    cb.set_title('Log-contact\ncorrelation',fontsize=10,pad=10)
    fig.suptitle('TDP-43: cross-temperature state matching',fontsize=18,y=.98)
    fig.text(.5,.025,'Boxed values: accepted one-to-one matches at correlation >= 0.90. Local indices shown one-based.',
             ha='center',fontsize=9)
    for ext in ['png','pdf']:fig.savefig(OUT/f'TDP43_matching_after_fix.{ext}',dpi=200)
    plt.close(fig)


def write_report(rows,contrasts):
    validation=json.loads((OUT/'contact_rebuild_validation.json').read_text())
    lines=['TDP-43 接触图索引修复：仅用于 R2C2 回复的独立分析',
           '', '修改范围：本目录。没有覆盖原状态映射、正文/SI 数据、绘图脚本或图文件；没有启动训练或 MD。',
           '高温已有正确的 148 残基索引，其接触图与模型种子对应原样保留。',
           '两套低温数据从完整距离特征与各模型的 soft assignments 重建，分别包含 20 个模型、200 万帧。', '']
    for r in validation:
        if r.get('rebuild_performed'):
            lines.append(f"{r['condition']}: 平均图复现最大误差 {r['old_average_error'][0]:.3g}; "
                         f"筛帧方式 {r['mode']}; 原种子对应全部复现; 修复后改变对应的模型数 {r['fixed']['changed_models']}/20。")
    lines+=['', '修复内容：',
            '正确残基对为 148 残基、序列间隔 >2。旧低温图误按 146 残基的全部上三角排列相同长度的特征向量。',
            '使用正确残基对重建后，沿用原同温度 constrained k-medoids/Hungarian 对齐程序及 gap=4 遮罩。',
            '跨温度仍按占比 >=0.10、log-contact correlation >=0.90 做一对一匹配，未为恢复特定结论而调阈值。',
            '这个阈值来自现有代码与原图对应关系的复现，不等于找回了全部历史运行环境日志。','',
            '主要结果：共同状态的 population-weighted transition-matrix distance（lag=2 ns）。',
            '所有比较都使用完整的未筛边转移概率；“原映射”不是正文原始绘图边筛选距离。']
    for variant,label in [('original','原映射'),('index_only','仅索引修复'),('fully_fixed','索引修复并重做种子对齐')]:
        f=previous.old.one(rows,variant=variant,cutoff=.90,family='free',uncertainty='sampling',metric='D_common_pooled_weights')
        d=previous.old.one(rows,variant=variant,cutoff=.90,family='condensate',uncertainty='sampling',metric='D_common_pooled_weights')
        delta=previous.old.one(contrasts,variant=variant,cutoff=.90,uncertainty='sampling',metric='D_common_pooled_weights')
        model=previous.old.one(contrasts,variant=variant,cutoff=.90,uncertainty='model_seeds',metric='D_common_pooled_weights')
        lines.append(f"{label}: free={f['estimate']:.6f}; condensate={d['estimate']:.6f}; "
                     f"cond-free={delta['estimate']:+.6f}; 采样95%区间=[{delta['ci_low']:+.6f},{delta['ci_high']:+.6f}]; "
                     f"模型种子95%区间=[{model['ci_low']:+.6f},{model['ci_high']:+.6f}]。")
    lines+=['','修复后的共同状态（以下 local state 从 1 开始，不是沿用正文所有全局编号）：']
    for family in s.FAMILIES:
        _,corr,accepted,_=match(family,'fully_fixed',.90)
        text='; '.join(f'TL local{lo+1} <-> TH local{hi+1}, r={corr[hi,lo]:.6f}' for lo,hi in accepted)
        lines.append(f'{family}: {text}')
    lines+=['','跨温度相关度阈值敏感性（占比筛选固定 0.10）：']
    for cutoff in [.80,.85,.90,.95]:
        delta=previous.old.one(contrasts,variant='fully_fixed',cutoff=cutoff,uncertainty='sampling',metric='D_common_pooled_weights')
        lines.append(f"阈值 {cutoff:.2f}: cond-free={delta['estimate']:+.6f}, 采样95%区间 "
                     f"[{delta['ci_low']:+.6f},{delta['ci_high']:+.6f}]")
    point=previous.old.one(contrasts,variant='fully_fixed',cutoff=.90,uncertainty='sampling',metric='D_common_pooled_weights')
    lines+=['','解释：',
            ('修复后的点估计仍是 condensate < free。' if point['estimate']<0 else '修复后的点估计是 condensate > free。'),
            '点估计方向不等于统计上确定的差异，应同时看独立列出的两类区间和匹配阈值敏感性。',
            '本分析回答共同状态之间的固定 lag 转移变化，不能概括全部非共同状态的物理行为。',
            '未匹配不等于物理不存在；不能把索引导致的匹配失败解释成状态出现/消失。',
            '仍是各条件独立训练后匹配，不是四条件共享同一个状态函数。',
            'Bootstrap 固定匹配关系，没有覆盖所有状态定义/匹配不确定性；没有把有限 lag 转移概率叫作连续时间速率。',
            '采样区间与模型种子区间各 4000 次重复，分开报告，不视为总误差区间。',
            'Free 以完整轨迹重采样；condensate 以绝对 slab 时间块重采样，主块长50 ns，另检查10/100 ns。',
            '', '文件：', 'TDP43_before_after.png/pdf：三种处理的距离及差值。',
            'TDP43_matching_after_fix.png/pdf：修复前后接触图相关性及被接受的匹配。',
            'distance_results.csv、condensate_minus_free.csv：全部距离、占比指标及区间。',
            'fixed_population_changes.csv、fixed_transition_entry_changes.csv：修复后的逐状态/逐转移变化。',
            'fixed_full_local_matrices.csv：完整五状态矩阵，不对缺失状态填零。',
            'source_manifest_before.json、original_figures_before.json：原始输入与图文件校验记录。',
            'SCOPE.txt：本次修改仅对独立 R2C2 回复分析生效，不回写正文或 SI。']
    (OUT/'RESULTS_zh.txt').write_text('\n'.join(lines)+'\n')


def main():
    rng=np.random.default_rng(SEED)
    data={}
    for family in s.FAMILIES:
        for role in s.ROLES:
            print('Resampling',family,role,flush=True)
            data[family,role]=sample_condition(family,role,rng)
    rows,contrasts=analyze(data)
    plot_comparison(rows,contrasts)
    plot_matching()
    write_report(rows,contrasts)
    provenance=dict(seed=SEED,nboot=NBOOT,lag_ns=2.,population_cutoff=.10,correlation_cutoff=.90,
        sensitivity_cutoffs=[.80,.85,.90,.95],primary_slab_block_ns=50,slab_block_sensitivity_ns=[10,100],
        mapping_bootstrapped=False,training_performed=False,traditional_MSM_used=False,
        note='Conditional on the specified matching procedure; not a four-condition shared state function')
    (OUT/'comparison_settings.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('DONE',OUT,flush=True)


if __name__=='__main__':main()
