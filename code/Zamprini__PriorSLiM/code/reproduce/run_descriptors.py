from pathlib import Path
import argparse
import os
import sys


def effects(root, output):
    import descriptor_effects as a
    from scipy.special import expit
    import numpy as np
    import pandas as pd
    source=root/'prior_compression_experiments/MOTIFAE_ONLY_COMPARISON/CV10_TEST/results'
    first=a.table(source/'fold_01/test_metadata.tsv')
    features, first=a.phenotype(first)
    rows=[]
    for fold in range(1,11):
        folder=source/f'fold_{fold:02d}'
        contract=a.t.read(folder/'evaluation_complete.json')
        for name, digest in contract['files'].items():
            a.t.require(a.t.digest(folder/name)==digest,'Prediction input changed')
        meta=a.table(folder/'test_metadata.tsv')
        a.t.require(np.array_equal(meta.window_id,first.window_id),'Test alignment')
        y=meta.label.to_numpy(int)
        for role in ['performance','efficiency']:
            zf,ft=a.load_prediction(folder/f'Full_{role}_test.npz',meta)
            for baseline in ['Explicit_only','MotifAE_only']:
                zb,bt=a.load_prediction(folder/f'{baseline}_test.npz',meta)
                b,f=expit(zb)>=bt,expit(zf)>=ft
                for name,label in [('positive',1),('negative',0)]:
                    case=(y==label)&(b!=y)&(f==y)
                    control=(y==label)&(b!=y)&(f!=y)
                    for feature,_ in a.EFIELDS:
                        values=features[feature].to_numpy();sd=float(values[y==label].std())
                        difference=float(values[case].mean()-values[control].mean()) if case.any() and control.any() else np.nan
                        rows.append(dict(fold=fold,role=role,baseline=baseline,calibration='own_val_F1',analysis=name,
                            feature=feature,case_n=int(case.sum()),control_n=int(control.sum()),
                            standardized_difference=difference/sd if sd>0 else np.nan,reference_sd=sd))
    frame=pd.DataFrame(rows);frame.to_csv(output/'e_effects.tsv',sep='\t',index=False)
    original=a.table(source/'summary/e_effects.tsv')
    keys=['fold','role','baseline','calibration','analysis','feature']
    merged=frame.merge(original,on=keys,suffixes=('_new','_old'),validate='one_to_one')
    a.t.require(len(merged)==len(frame)==len(original),'Descriptor regression row count')
    for field in ['standardized_difference','reference_sd','case_n','control_n']:
        np.testing.assert_allclose(merged[field+'_new'],merged[field+'_old'].astype(float),atol=1e-12,rtol=1e-10,equal_nan=True)
    print(f'PASS: {len(frame)} descriptor-effect rows match recorded results')


def composition(root,output,extract):
    import composition as a
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    a.torch.set_num_threads(1)
    source=root/'prior_compression_experiments/MOTIFAE_EVIDENCE_DISCOVERY/results'
    if extract:
        meta,_=a.metadata()
        dev=pd.concat([a.prediction_rows(meta,fold,'efficiency','development')[0] for fold in range(1,11)],ignore_index=True)
        a.write(dev,output/'development_window_transitions.tsv')
        mat,position,selected=a.sparse_contributions(dev,output)
    else:
        dev=a.table(source/'development_window_transitions.tsv')
        mat=sp.load_npz(source/'development_contributions_by_feature.npz')
        selected=np.load(source/'development_contributions_by_position.npz')['selected_features_by_fold']
    for column in ['discovery_eligible','canonical29']:
        dev[column]=dev[column].astype(str).eq('True')
    a.require(set(dev.cohort)=={'development'},'Expected development data')
    existing,new,definitions=a.biological_features(dev)
    signatures=a.negative_signatures(dev,mat,selected)
    discovery=a.association_tests(dev,existing,new,mat,selected,signatures,'discovery',cohort='motif_matched')
    a.write(discovery,output/'negative_discovery_associations.tsv')
    candidates=discovery[discovery.descriptor.str.startswith('aa_')&(discovery.q<.05)&(abs(discovery.partial_r)>=.15)].sort_values('q')
    locked=candidates.groupby('feature_id',sort=False).head(2).head(12)
    pairs=pd.DataFrame({'feature_id':locked.feature_id.unique()})
    rep=a.association_tests(dev,existing,new,mat,selected,pairs,'replication',locked,cohort='motif_matched')
    rep=rep.merge(locked[['feature_id','descriptor','partial_r','q']].rename(columns={'partial_r':'discovery_partial_r','q':'discovery_q'}),on=['feature_id','descriptor'],validate='one_to_one')
    rep['same_direction']=np.sign(rep.partial_r)==np.sign(rep.discovery_partial_r)
    rep['replicated']=(rep.q<.05)&rep.same_direction
    a.write(rep,output/'fine_composition_replication.tsv')
    expected=a.table(source/'fine_composition_replication.tsv')
    merged=rep.merge(expected,on=['feature_id','descriptor'],suffixes=('_new','_old'),validate='one_to_one')
    a.require(len(merged)==len(rep)==len(expected),'Composition regression row count')
    for field in ['partial_r','p','q','discovery_partial_r','discovery_q']:
        np.testing.assert_allclose(merged[field+'_new'],merged[field+'_old'].astype(float),atol=1e-9,rtol=1e-6)
    print(f'PASS: {len(discovery)} discovery tests and {len(rep)} follow-up tests; follow-up matches recorded results')


def main():
    p=argparse.ArgumentParser(description='Recompute corrected-window descriptor effects and composition associations.')
    p.add_argument('phase',choices=['effects','composition','composition-extract'])
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();root=args.data_root.resolve();output=args.output_dir.resolve()
    if output.is_relative_to(root):p.error('Output must be outside input data')
    if output.exists():p.error('Choose a new output directory')
    os.environ['PRIORSLIM_DATA_ROOT']=str(root)
    here=Path(__file__).resolve().parent
    sys.path[:0]=[str(here.parent/'src'),str(here/'ablation_helpers')]
    output.mkdir(parents=True)
    if args.phase=='effects':effects(root,output)
    else:composition(root,output,args.phase=='composition-extract')

if __name__=='__main__':main()
