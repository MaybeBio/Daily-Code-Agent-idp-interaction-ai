#!/usr/bin/env python3
"""Concise bilingual analysis summaries, not edits to the response letter."""
import csv
import json
import numpy as np
from extract import HERE,CONFIG,dump
from analyze import table


def read(name):
    with (HERE/name).open() as f:return list(csv.DictReader(f))


def value(row):
    if row['ratio']=='':return 'not resolved'
    result=f"{float(row['ratio']):.3f}"
    if row['status']=='finite':result+=f" [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}]"
    return result


def main():
    mobility=read('mobility_summary.csv');rr=read('temperature_ratios.csv');counts=read('residence_counts.csv')
    occ=read('vamp_sampled_residence.csv');prov=json.loads((HERE/'training_provenance.json').read_text())
    summary=[];chainids=[]
    for case in CONFIG:
        primary=next(r for r in mobility if r['case']==case and r['lag_ns']=='10' and r['coordinates']=='xy' and
            r['selected']=='original_window_matched')
        spatial=next(r for r in mobility if r['case']==case and r['lag_ns']=='10' and r['coordinates']=='xy' and
            r['selected']=='original_ids_bulk_origins')
        c=[r for r in counts if r['case']==case];o=[r for r in occ if r['case']==case]
        pr=next(r for r in prov if r['case']==case)
        denom=sum(int(r['observed_frames']) for r in o)
        frac=sum(int(r['observed_frames'])*float(r['source_rule_bulk_fraction']) for r in o)/denom
        uniform=sum(int(r['observed_frames'])*float(r['uniform_rule_bulk_fraction']) for r in o)/denom
        summary.append(dict(case=case,original_inter_chains=pr['inter_used_chains'],
            traced_VAMP_trajectories=pr['coordinate_verified_trajectories'],
            additional_untraced_VAMP_trajectories=len(pr['untraced_training_trajectories']),
            matched_window_mobility_ratio=float(primary['ratio']),mobility_ci_low=float(primary['ci_low']),
            mobility_ci_high=float(primary['ci_high']),matched_window_blocks=int(primary['n_blocks']),
            bulk_matched_mobility_ratio=float(spatial['ratio']),
            min_long_chains=min(int(r['n_500']) for r in c),max_long_chains=max(int(r['n_500']) for r in c),
            min_short_chains=min(int(r['n_100']) for r in c),max_short_chains=max(int(r['n_100']) for r in c),
            sampled_VAMP_chains=len({r['chain_id'] for r in o}),sampled_VAMP_chain_frames=denom,
            sampled_VAMP_source_rule_bulk_fraction=frac,sampled_VAMP_uniform_rule_bulk_fraction=uniform))
        for path in sorted((HERE/'raw').glob(f'{case}_block*.npz')):
            with np.load(path) as z:
                flags=z['bulk'];lo=flags.all(axis=0);sh=flags.reshape(5,500,200).all(axis=1)
                chainids.append(dict(case=case,block=int(path.stem.rsplit('block',1)[1]),
                    first_time_ns=float(z['time_ns'][0]),last_time_ns=float(z['time_ns'][-1]),
                    long_chain_ids=np.flatnonzero(lo).tolist(),short_chain_ids=[np.flatnonzero(x).tolist() for x in sh]))
    table('selection_summary.csv',summary);dump(HERE/'new_residence_chain_ids.json',chainids)
    zh=['R2C5 分析结果（尚未合并正式回复信、正文或 SI）','',
        '范围：六个凝聚体条件；每个条件四个分散的 500 ns 窗口，共 2 us，全部 200 条链。',
        '这是分布于生产段的抽样对照，不是整个生产轨迹逐帧重算。原始时间间隔已核对。','',
        '1. 原入选链相对全部链的运动性',
        '主比较：10 ns 横向质量质心 MSD；严格匹配原选链窗口中的时间原点。',
        '比值小于 1 表示入选链比全部 200 条链的平均运动更慢。括号内为时间块重采样区间。']
    en=['R2C5 ANALYSIS RESULTS (not integrated into the formal revision)','',
        'Six condensate conditions; four distributed 500 ns windows per condition, 2 us total per condition, all 200 chains.',
        'This is a distributed production-trajectory subsample, not an exhaustive whole-production reanalysis.','',
        '1. Mobility: selected/all lateral mass-COM MSD at 10 ns, matched to original source-window origins.']
    for r in summary:
        line=f"{r['case']}: {r['matched_window_mobility_ratio']:.3f} [{r['mobility_ci_low']:.3f}, {r['mobility_ci_high']:.3f}]; bulk-matched ratio {r['bulk_matched_mobility_ratio']:.3f}; source-window blocks {r['matched_window_blocks']}"
        zh.append(line);en.append(line)
    zh+=['','不能把“选中链相对全部链较慢”直接等同于筛选制造了接触松弛趋势。全部链也包含界面/稀相链；',
         '因此另列空间匹配的中心区比较。是否影响接触动力学，要看下面真正缩短驻留要求的对照。','',
         '2. 原入选片段与 100 ns 驻留筛选的直接比较',
         '原入选链必须完整覆盖相应 100 ns 子段；新筛选组在完全相同的子段和父时间块计算。',
         '空间标准按现存条件脚本保持一致。两组均用 100 ns 子段计算 ACF/伙伴保留，避免片段长度混杂。',
         '以下为高温/低温阈值时间比值：接触 ACF 阈值 0.2，伙伴保留阈值 0.5。']
    en+=['','The all-chain reference includes interface/dilute chains. A spatially matched bulk-origin comparison is therefore also provided.',
         'A mobility difference alone does not establish an artifact in contact relaxation.','',
         '2. Historical selected windows versus strict 100 ns residence: both groups use identical 100 ns analysis segments and supported parent blocks.',
         'The historical group is defined by the original metadata and retained chain IDs, not by an assumed perfect reconstruction of the historical bulk fraction.',
         'Heating ratios of first ensemble-curve crossing times: contact ACF at 0.2; partner retention at 0.5.']
    for system in ['FUS','A1','TDP43']:
        for metric in ['inter','partner']:
            rows=[r for r in rr if r['system']==system and r['metric']==metric and r['quantity']=='tau' and
                r['group'] in ['original_selection','short_matched_original']]
            labels={'original_selection':'original','short_matched_original':'100 ns'}
            line=f"{system} {metric}: "+'; '.join(labels[r['group']]+': '+value(r) for r in rows)
            zh.append(line);en.append(line)
    short=[r for r in rr if r['group']=='short_matched_original' and r['quantity']=='tau']
    directional=len(short)==6 and all(r['ratio']!='' and (float(r['ratio'])>1 if r['metric']=='inter' else float(r['ratio'])<1) for r in short)
    if directional:
        zh+=['','三个体系在 100 ns 驻留规则下的点估计仍为：接触数阈值松弛变慢、伙伴更快更换。',
             '因此，在本次抽样窗口和已声明的空间筛选规则下，这两个方向不只存在于原选链集合中。',
             '这不是“完全没有选择偏差”或“所有误差范围均不跨 1”的声明；区间必须逐体系、逐指标阅读。']
        en+=['','Under the 100 ns rule, point estimates retain slower contact-number threshold relaxation and faster partner exchange upon heating in all three sequences.',
             'Within these sampled windows and stated spatial rules, those directions are not exclusive to the historically selected cohort.',
             'This does not establish absence of selection bias or that every interval excludes unity.']
    else:
        zh+=['','各体系结果需分别解读，不能概括为全部体系均通过同方向的驻留稳健性检验。']
        en+=['','Interpret the sequences separately; do not claim universal directional robustness.']
    zh+=['','另做严格 500 ns 与严格 100 ns 重选对照，数值保存在 CSV 的 residence_500 / residence_100 中。',
         '若某窗口严格 500 ns 合格链为零，只报告缺失，不能冒充已完成该窗口的配对比较。',
         '全部四个窗口的 100 ns 结果另存 residence_100_all_blocks。',
         '严格重选合格链数量（每个 500 ns 父窗口 / 每个 100 ns 子窗口）：']
    en+=['','A separate strict 500/100 ns re-selection control is retained in CSV groups residence_500/residence_100.',
         'A window with no strict 500 ns resident is missing from that paired control, never a fabricated zero curve.',
         'The unrestricted four-block 100 ns result is in residence_100_all_blocks.',
         'Strictly eligible chain ranges per 500 ns parent / per 100 ns subwindow:']
    for r in summary:
        line=f"{r['case']}: {r['min_long_chains']}--{r['max_long_chains']} / {r['min_short_chains']}--{r['max_short_chains']}"
        zh.append(line);en.append(line)
    zh+=['','3. VAMPNet 选链及界面问题',
        '重新核对了每个条件 40 条主要训练轨迹，共 240 条；首、中、末共 720 帧与原始父链周期坐标完全一致。',
        '这些 VAMPNet 轨迹与链间分析来自同一套父选链窗口，不是直接取全部 200 条链。29--40 的差异还涉及分析段范围过滤。',
        '现存提取脚本：FUS 高温用 5 nm、10--90% 中心区；A1 高温用 4.5 nm、20--80%；其余用 4 nm、20--80%。',
        'TDP-43 的准备脚本默认允许 80% 驻留；原运行参数覆盖未归档，不能把这个默认值当作已证实的历史命令。',
        'A1 高温另有 chain139 的 7794 帧训练输入缺绝对来源窗口，不能认证其 bulk/interface 来源；本次未重新训练模型。',
        '因此不能写“所有训练片段的每一帧都严格排除了界面”。可写其来源为 bulk-oriented 选链，同时说明实际规则和例外。',
        '本次 sampled residence 表记录的是操作性中心区分类，不是单独重建的液气物理界面标签。','',
        '原训练链在此次可观测窗口内的中心区占比（现存脚本标准 / 统一 4 nm、20--80% 标准）：']
    en+=['','3. Training-chain provenance and the interface question',
        'All 240 principal training DCDs (40 per condition) were checked at their first, middle and last frames: 720 exact periodic-coordinate matches to raw parent chains.',
        'They share parent selected windows with the intermolecular analysis; they are not an unrestricted sample of all 200 chains.',
        'Extant spatial rules differ: FUS high T uses 5 nm and 10--90%; A1 high T uses 4.5 nm and 20--80%; the remaining cases use 4 nm and 20--80%.',
        'The TDP-43 preparation wrapper defaults to 80% residence; historical environment/CLI overrides are not archived, so the default is not certified as the actual invocation.',
        'One additional A1 high-T training input, chain139 with 7,794 frames, lacks its absolute source-window metadata. Its bulk/interface provenance remains uncertified; no model was retrained here.',
        'A bulk-oriented source selection is supported. Exclusion of every interfacial training frame is not supported.',
        'Percentile/cluster exits below are operational central-region exits, not a separately reconstructed physical liquid-vapor interface label.','',
        'Sampled training-window bulk fraction: extant spatial rule / uniform 4 nm, 20--80% rule:']
    for r in summary:
        line=f"{r['case']}: {r['sampled_VAMP_source_rule_bulk_fraction']:.4f} / {r['sampled_VAMP_uniform_rule_bulk_fraction']:.4f}; sampled chains {r['sampled_VAMP_chains']}"
        zh.append(line);en.append(line)
    zh+=['','统计边界：4,000 次重采样以 500 ns 整个时间块为单位，同一 slab 的链不当作独立模拟重复。',
        '每条件仅四个抽样块；部分原选链窗口比较使用更少的重叠块。区间是条件性时间块不确定度，不是独立重复证明。',
        '正式回复、正文和 SI 均未改动。可先看 All_C5_analysis_figures.pdf；定义及复现入口见 METHODS.txt。']
    en+=['','Uncertainty: 4,000 resamples of whole 500 ns temporal blocks, not independent chain replicas.',
        'Only four sampled blocks per condition; some source-window comparisons have fewer overlapping blocks. Intervals are conditional temporal-block uncertainty, not independent-run reproducibility.',
        'Formal response letter, manuscript and SI are unchanged. See All_C5_analysis_figures.pdf and METHODS.txt.']
    (HERE/'RESULTS_ZH.txt').write_text('\n'.join(zh)+'\n')
    (HERE/'RESULTS_EN.txt').write_text('\n'.join(en)+'\n')


if __name__=='__main__':main()
