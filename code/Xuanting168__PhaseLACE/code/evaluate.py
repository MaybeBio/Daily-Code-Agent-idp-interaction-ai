from __future__ import annotations

import argparse
import csv
import glob
import json
from pathlib import Path
from statistics import mean

import numpy as np
from sklearn.metrics import average_precision_score,auc,brier_score_loss,precision_recall_curve,roc_auc_score


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    if y.ndim!=1 or y.shape!=p.shape or set(np.unique(y))!={0,1}:
        raise ValueError('Evaluation requires matching one-dimensional labels with both classes')
    if not np.isfinite(p).all() or np.any((p<0)|(p>1)):raise ValueError('Invalid probabilities')
    precision,recall,_=precision_recall_curve(y,p)
    return {'auroc':float(roc_auc_score(y,p)),'aupr':float(auc(recall,precision)),
            'ap':float(average_precision_score(y,p)),'brier':float(brier_score_loss(y,p))}


def fold_metrics(y: np.ndarray, p: np.ndarray, fold: int) -> dict:
    result=metrics(y,p);pos=np.flatnonzero(y==1);neg=np.flatnonzero(y==0)
    size=int(.8*min(len(pos),len(neg)))
    if size<1:raise ValueError('Not enough samples for balanced80 evaluation')
    rng=np.random.RandomState(2026091600+fold);balanced={k:0. for k in ['auroc','aupr','ap']}
    for _ in range(100):
        negative=rng.choice(neg,size=size,replace=False)
        positive=rng.choice(pos,size=size,replace=False)
        ids=np.r_[negative,positive];values=metrics(y[ids],p[ids])
        for k in balanced:balanced[k]+=values[k]/100
    result.update({'balanced80_'+k:v for k,v in balanced.items()})
    return result


def evaluate(paths: list[Path], members: int, index_path: Path | None = None) -> dict:
    if not paths or len({p.resolve() for p in paths})!=len(paths):raise ValueError('Missing or duplicate prediction files')
    index={}
    if index_path:
        with index_path.open(newline='',encoding='utf-8-sig') as h:
            original=list(csv.DictReader(h));index={int(r['row']):r for r in original}
        if len(index)!=len(original):raise ValueError('Duplicate row in index')
    combined={};per_fold_members={};per_fold_keys={}
    for path in paths:
        with path.open(newline='',encoding='utf-8-sig') as h:rows=list(csv.DictReader(h))
        if not rows:raise ValueError('Empty prediction file')
        seen=set();this_fold_keys={}
        for row in rows:
            i=int(row['row']);source=index.get(i)
            if index and source is None:raise ValueError('Prediction row absent from index')
            if source and any(row.get(k)!=source.get(k) for k in ['label','sequence_sha256','accession','original_row']):
                raise ValueError('Prediction identity does not match index')
            fold=int(row['fold'] if 'fold' in row else source['fold']) if ('fold' in row or source) else None
            if fold is None or fold not in range(5):raise ValueError('Supply a fold column or --index with frozen fold membership')
            if source and fold!=int(source['fold']):raise ValueError('Fold differs from index')
            key=(fold,i)
            if key in seen:raise ValueError('Repeated sample within one prediction file')
            seen.add(key);this_fold_keys.setdefault(fold,set()).add(i)
            p=float(row['probability']);label=int(row['label'])
            if label not in (0,1) or not np.isfinite(p) or not 0<=p<=1:raise ValueError('Invalid label/probability')
            identity=tuple(row.get(k,'') for k in ['label','sequence_sha256','accession','original_row'])
            if key not in combined:combined[key]={'identity':identity,'label':label,'p':[]}
            if combined[key]['identity']!=identity:raise ValueError('Member identities differ')
            combined[key]['p'].append(p)
        for fold,keys in this_fold_keys.items():
            if fold in per_fold_keys and per_fold_keys[fold]!=keys:raise ValueError('Member sample coverage differs')
            per_fold_keys[fold]=keys;per_fold_members[fold]=per_fold_members.get(fold,0)+1
    if any(n!=members for n in per_fold_members.values()) or any(len(r['p'])!=members for r in combined.values()):
        raise ValueError('Every sample/fold must have exactly the requested number of members')
    values=[]
    for fold in sorted(per_fold_keys):
        rows=[combined[(fold,i)] for i in sorted(per_fold_keys[fold])]
        y=np.array([r['label'] for r in rows]);p=np.array([mean(r['p']) for r in rows])
        values.append({'fold':fold,'n':len(rows),'positives':int(y.sum()),**fold_metrics(y,p,fold)})
    keys=[k for k in values[0] if k not in ['fold','n','positives']]
    return {'members':members,'outer_folds':len(values),'folds':values,'mean_fold_metrics':{k:mean(r[k] for r in values) for k in keys},
            'protocol':'Equal member probability average; metrics within fold, then equal fold average; balanced80 uses 100 paired draws.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Evaluate final-model prediction CSVs; wildcard patterns are expanded internally.')
    parser.add_argument('--predictions',nargs='+',required=True)
    parser.add_argument('--members',type=int,choices=[1,3,6],default=3)
    parser.add_argument('--index',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();paths=[]
    for pattern in args.predictions:
        found=sorted(glob.glob(pattern))
        if not found:raise FileNotFoundError(pattern)
        paths.extend(Path(p) for p in found)
    result=evaluate(paths,args.members,args.index)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result['mean_fold_metrics']))
