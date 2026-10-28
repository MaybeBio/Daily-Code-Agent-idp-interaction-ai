from __future__ import annotations
import numpy as np
from scipy.stats import spearmanr

def region_window(scores: np.ndarray, threshold: float, minimum: int) -> np.ndarray:
    """Reproduce author's compute_window marking and gap filling in linear time."""
    active = np.asarray(scores) > threshold
    n = len(active)
    marks = np.zeros(n, dtype=np.int8)
    bounds = np.flatnonzero(np.diff(np.r_[False, active, False].astype(int)))
    for start, end in zip(bounds[::2], bounds[1::2]):
        if end - start >= minimum:
            marks[start:min(end + 1, n)] = 1
    gaps = np.flatnonzero(np.diff(np.r_[False, marks == 0, False].astype(int)))
    for start, end in zip(gaps[::2], gaps[1::2]):
        if end - start < minimum:
            marks[start:end] = 1
    return marks

def metric(y: np.ndarray, p: np.ndarray, intervals: list[tuple], predicted: np.ndarray) -> dict:
    informative = np.unique(y).size == 2 and int((y == 0).sum()) >= 20
    rho = float(spearmanr(y, p).statistic) if informative and np.unique(p).size > 1 else None
    overlap = int(np.sum((y == 1) & (predicted == 1)))
    hit = sum((bool(predicted[a:b].any()) for a, b in intervals))
    return dict(spearman=rho, informative=informative, annotated_residues=int(y.sum()), predicted_residues=int(predicted.sum()), overlap=overlap, regions=len(intervals), hit_regions=hit)

import argparse
import json
from pathlib import Path

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Evaluate one local score vector against 1-based closed annotated intervals.')
    parser.add_argument('--scores',type=Path,required=True,help='One-dimensional NPY of 33-aa local probabilities')
    parser.add_argument('--annotations',type=Path,required=True,help='JSON list of [start,end] pairs, 1-based closed')
    parser.add_argument('--coordinate-mode',choices=['author','standard'],default='author')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();p=np.load(args.scores,allow_pickle=False)
    intervals=json.loads(args.annotations.read_text());y=np.zeros(len(p),dtype=np.int8);slices=[]
    if p.ndim!=1 or not len(p) or not np.isfinite(p).all() or np.any((p<0)|(p>1)):
        raise ValueError('Scores must be a finite probability vector')
    if not intervals:raise ValueError('At least one annotated interval is required')
    for a,b in intervals:
        if not isinstance(a,int) or not isinstance(b,int) or not 1<=a<=b<=len(p):raise ValueError('Invalid 1-based closed interval')
        start=a if args.coordinate_mode=='author' else a-1;y[start:b]=1;slices.append((start,b))
    selected=region_window(region_window(p,.5,5),.9,20)
    result={'coordinate_mode':args.coordinate_mode,**metric(y,p,slices,selected)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result))
