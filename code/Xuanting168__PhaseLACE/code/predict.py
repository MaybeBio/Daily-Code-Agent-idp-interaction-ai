from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from safetensors import safe_open
from safetensors.torch import load_file

from model import BaselineHead,PhaseLACE,pool_windows


def load_model(path: Path, device: str = 'cpu') -> tuple[torch.nn.Module,dict]:
    if path.suffix=='.safetensors':
        state=load_file(str(path),device='cpu')
        with safe_open(str(path),framework='pt',device='cpu') as handle:metadata=handle.metadata() or {}
    else:
        payload=torch.load(path,map_location='cpu',weights_only=True)
        if not isinstance(payload,dict):raise ValueError('Expected a tensor state dict or inference.pt')
        state=payload.get('model',payload)
        metadata={k:v for k,v in payload.items() if k!='model'} if 'model' in payload else {}
    if not state or not all(isinstance(k,str) and isinstance(v,torch.Tensor) for k,v in state.items()):
        raise ValueError('Use the tensor-only inference.pt exported by train.py')
    residual=any(k.startswith('baseline.') for k in state)
    model=PhaseLACE(BaselineHead(),float(metadata.get('residual_scale',.3))) if residual else BaselineHead()
    model.load_state_dict(state,strict=True)
    return model.to(device).eval(),metadata


def predict(features: np.ndarray, checkpoints: list[Path], baseline_checkpoints: list[Path] | None = None,
            device: str = 'cpu') -> dict:
    if features.ndim!=2 or features.shape[1]!=650 or not len(features) or not np.isfinite(features).all():
        raise ValueError('Expected finite nonempty [L,650] ESM320+ALBATROSS330 features')
    baseline_checkpoints=baseline_checkpoints or []
    if len(checkpoints)!=3 or len(baseline_checkpoints) not in (0,3):
        raise ValueError('PhaseLACE uses 3 residual members; Fusion additionally uses 3 base members')
    paths=checkpoints+baseline_checkpoints
    if len({p.resolve() for p in paths})!=len(paths):raise ValueError('Duplicate checkpoints')
    x=torch.tensor(features,dtype=torch.float32,device=device);probabilities=[];curves=[];folds=set()
    with torch.inference_mode():
        for i,path in enumerate(paths):
            model,metadata=load_model(path,device)
            if 'fold' in metadata:folds.add(int(metadata['fold']))
            if isinstance(model,PhaseLACE):
                if i>=3:raise ValueError('Fusion base group contains a residual member')
                parts=model.forward_components(x);p=parts['probability']
                logits=parts['base_position_logits']+model.residual_scale*parts['delta_position_logits']
            else:
                if i<3:raise ValueError('Main group requires PhaseLACE residual members')
                pooled=pool_windows(x);p=model(pooled);logits=model.position_logits(pooled)
            probabilities.append(float(p.cpu()));curves.append(logits[0].sigmoid().cpu().numpy())
    if len(folds)>1:raise ValueError('Members belong to different outer-fold models')
    return {'model':'PhaseLACE-Fusion' if baseline_checkpoints else 'PhaseLACE',
            'fold':next(iter(folds)) if folds else None,'probability':sum(probabilities)/len(probabilities),
            'member_probabilities':probabilities,'residue_33aa_mean':np.mean(curves,axis=0).tolist()}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Predict using external PhaseLACE checkpoints; no weights are bundled.')
    parser.add_argument('--features',type=Path,required=True)
    parser.add_argument('--checkpoints',type=Path,nargs=3,required=True)
    parser.add_argument('--baseline-checkpoints',type=Path,nargs=3)
    parser.add_argument('--device',default='cpu')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();torch.set_num_threads(4)
    result=predict(np.load(args.features,allow_pickle=False),args.checkpoints,args.baseline_checkpoints,args.device)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
