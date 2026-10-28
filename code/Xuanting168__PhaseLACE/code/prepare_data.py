from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import torch
from sklearn.model_selection import KFold

from model import pool_windows
from utils import sha256_file


def prepare(table: Path, features_dir: Path, output: Path, device: str = 'cpu') -> None:
    with table.open(newline='',encoding='utf-8-sig') as handle:
        source=list(csv.DictReader(handle))
    required={'feature_file','accession','label','category','sequence_sha256'}
    if not source or not required<=set(source[0]):
        raise ValueError('Required CSV fields: '+', '.join(sorted(required)))
    folds=np.full(len(source),-1,dtype=int)
    if 'fold' not in source[0]:
        for fold,(_,ids) in enumerate(KFold(5,shuffle=True,random_state=711016).split(source)):folds[ids]=fold
    rows=[];paths=[];offset=0
    for i,r in enumerate(source):
        path=features_dir/r['feature_file'];x=np.load(path,mmap_mode='r',allow_pickle=False)
        if x.ndim!=2 or x.shape[1]!=650 or not len(x) or x.dtype!=np.float32:
            raise ValueError(f'Expected nonempty float32 [L,650]: {path.name}')
        label=int(r['label']);fold=int(r['fold']) if 'fold' in r else int(folds[i])
        if label not in (0,1) or fold not in range(5):raise ValueError('Invalid label or fold')
        if len(r['sequence_sha256'])!=64 or any(c not in '0123456789abcdef' for c in r['sequence_sha256']):
            raise ValueError('sequence_sha256 must be a lowercase SHA256 digest of the input sequence')
        if 'row' in r and int(r['row'])!=i:raise ValueError('Input row order differs from explicit row identifiers')
        rows.append({'row':i,'original_row':r.get('original_row',str(i)),'accession':r['accession'],
                     'category':r['category'],'label':label,'length':len(x),'offset':offset,'fold':fold,
                     'sequence_sha256':r['sequence_sha256']})
        offset+=len(x);paths.append(path)
        del x
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError('Use a new output directory')
    if shutil.disk_usage(output.parent).free < offset*650*4*4+512*1024**2:
        raise RuntimeError('Insufficient disk for raw and pooled caches')
    output.mkdir();raw_dir=output/'raw';pooled_dir=output/'pooled';raw_dir.mkdir();pooled_dir.mkdir()
    raw_path=raw_dir/'raw_fp32.npy';pooled_path=pooled_dir/'pooled_fp32.npy'
    raw=np.lib.format.open_memmap(raw_path,mode='w+',dtype=np.float32,shape=(offset,650))
    pooled=np.lib.format.open_memmap(pooled_path,mode='w+',dtype=np.float32,shape=(offset,3,650))
    torch.set_num_threads(4)
    for row,path in zip(rows,paths):
        x=np.load(path,allow_pickle=False)
        if not np.isfinite(x).all():raise ValueError('Non-finite features: '+path.name)
        start=row['offset'];end=start+row['length'];raw[start:end]=x
        with torch.inference_mode():pooled[start:end]=pool_windows(torch.as_tensor(x,device=device)).permute(1,0,2).cpu().numpy()
    raw.flush();pooled.flush();del raw,pooled
    for directory,path,shape in [(raw_dir,raw_path,[offset,650]),(pooled_dir,pooled_path,[offset,3,650])]:
        index=directory/'train_index.csv'
        with index.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        manifest={'cache_shape':shape,'cache_dtype':'float32','cache_bytes':path.stat().st_size,
                  'cache_sha256':sha256_file(path),'train_index_sha256':sha256_file(index),
                  'source_table_sha256':sha256_file(table)}
        (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Build raw and multiscale caches from externally supplied residue embeddings.')
    parser.add_argument('--table',type=Path,required=True)
    parser.add_argument('--features-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--device',default='cpu')
    args=parser.parse_args();prepare(args.table,args.features_dir,args.output,args.device)
