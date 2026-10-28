from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Literal, Any, Dict, Iterable
import yaml
import pandas as pd
import numpy as np

from idr_diverge.distances.compute_ndist import ComputeNDistanceDict, _load_embeddings, _load_ortholog_embeddings, _parse_ids,save_pickle, load_fam_map, calc_FND
from idr_diverge.utils.helpers import read_pickle, resolve_config_paths

#Command-line interface for computing embedding neighbour distances/divergence and family neighbour divergence (FND).
#Handles argument parsing and loading config file, and makes calls to embedding scripts in idr_diverge/compute_ndist
#This script was developed with assistance from ChatGPT (OpenAI) and has been reviewed and modified by Denise Le.

@dataclass(frozen=True)
class InputsConfig:
    dist_matrix: Optional[Path]
    background_emb: Optional[Path]
    ortholog_emb: Optional[Path]
    #homolog_annotation: Optional[Path]
    fam_map: Optional[Path]

@dataclass(frozen=True)
class OutputsConfig:
    out_distance: Path #directory path
    out_random: Path #file path
    #overwrite: bool = False

@dataclass(frozen=True)
class ComputeConfig:
    #metric: Literal["cosine", "euclidean", "manhattan"] = "cosine"
    chunk_size: int = 5000
    #n_jobs: int = 1
    job_id: int = 0

@dataclass(frozen=True)
class NeighbourDivergenceConfig:
    dist_type: Literal["between", "within"]
    segment: str = "IDR"
    return_mode: Literal["aggregate", "per_gene"] = "aggregate"

@dataclass(frozen=True)
class RandomNeighbourDivergenceConfig:
    random_seed: Optional[int]# = 42
    enabled: bool = False
    dist_type: Literal["between", "within"] = 'within'
    segment: str = "IDR"
    sample_size: int = 100
    
    return_mode: Literal["aggregate", "per_gene"] = "aggregate"

@dataclass(frozen=True)
class OrthologFilterConfig:
    enabled: bool = True
    homolog_annotation: Optional[Path] = None
    ortholog_types: Literal["ortholog_one2one","ortholog_one2many","ortholog_one2one_one2many"] = "ortholog_one2one_one2many"
    min_species_coverage: float = 0.5
    min_species: int = 10

@dataclass(frozen=True)
class FNDconfig:
    enabled: bool = False
    fam_distance_matrix: Optional[Path] = None
    bg_distance_matrix: Optional[Path] = None
    output_results: Path = Path("fnd_results.csv")

@dataclass(frozen=True)
class PipelineConfig:
    inputs: InputsConfig
    outputs: OutputsConfig
    compute: ComputeConfig
    neighbour_divergence: NeighbourDivergenceConfig
    random_neighbour_divergence: RandomNeighbourDivergenceConfig
    ortholog_filter: OrthologFilterConfig
    FND: FNDconfig
    raw: Dict[str, Any]  # keep raw dict for debugging/forward compatibility


def _to_path(x) -> Optional[Path]:
    if x is None:
        return None
    x = str(x).strip()
    return None if x.lower() in {"null", "none", ""} else Path(x)

def load_pipeline_config(path: str | Path) -> PipelineConfig:
    path = Path(path)
    with path.open("r") as f:
        cfg = yaml.safe_load(f)
    
    cfg = resolve_config_paths(config=cfg, config_path = path,
                               path_keys = {'background_emb', 'ortholog_emb', 'fam_map',
                                            'homolog_annotation', 'out_distance', 'out_random',
                                            'fam_distance_matrix', 'bg_distance_matrix', 'output_results'}
    )

    # inputs
    inputs = InputsConfig(
        dist_matrix=_to_path(cfg["inputs"].get("dist_matrix")),
        background_emb=_to_path(cfg["inputs"].get("background_emb")),
        ortholog_emb=_to_path(cfg["inputs"].get("ortholog_emb")),
        #homolog_annotation=_to_path(cfg["inputs"].get("homolog_annotation")),
        fam_map=_to_path(cfg["inputs"].get("fam_map")),
    )

    outputs = OutputsConfig(
        out_distance=Path(cfg["neighbour_divergence"]["out_distance"]),
        out_random = Path(cfg["random_neighbour_divergence"]["out_random"]),
        #overwrite=bool(cfg["outputs"].get("overwrite", False)),
    )

    compute = ComputeConfig(
        #metric=cfg["compute"].get("metric", "cosine"),
        chunk_size=int(cfg["compute"].get("chunk_size", 5000)),
        #n_jobs=int(cfg["compute"].get("n_jobs", 1)),
        job_id=int(cfg["compute"].get("job_id", 0))
    )

    neighbour_div = NeighbourDivergenceConfig(
        dist_type=cfg["neighbour_divergence"].get("dist_type", "between"),
        segment=cfg["neighbour_divergence"].get("segment", "IDR"),
        return_mode=cfg["neighbour_divergence"].get("return_mode", "aggregate"),
    )

    random_neighbour_div = RandomNeighbourDivergenceConfig(
        enabled=bool(cfg["random_neighbour_divergence"].get("enabled", False)),
        dist_type=cfg["neighbour_divergence"].get("dist_type", "within"),
        segment=cfg["neighbour_divergence"].get("segment", "IDR"),
        sample_size=int(cfg["random_neighbour_divergence"].get("sample_size", 100)),
        random_seed=(int(cfg["random_neighbour_divergence"]["random_seed"])
    if cfg["random_neighbour_divergence"].get("random_seed") is not None
    else None),
        return_mode=cfg["neighbour_divergence"].get("return_mode", "aggregate")
       
    )

    ortholog_filter = OrthologFilterConfig(
        enabled=bool(cfg["ortholog_filter"].get("enabled", True)),
        homolog_annotation=_to_path(cfg["ortholog_filter"].get("homolog_annotation")),
        ortholog_types=str(cfg["ortholog_filter"].get("ortholog_types", "ortholog_one2one_one2many")),
        min_species_coverage=float(cfg["ortholog_filter"].get("min_species_coverage", 0.5)),
        min_species=int(cfg["ortholog_filter"].get("min_species", 10)),
    )
    
    FND = FNDconfig( #TODO
        enabled=bool(cfg["calculate_FND"].get("enabled", False)),
        fam_distance_matrix=_to_path(cfg["calculate_FND"].get("fam_distance_matrix")),
        bg_distance_matrix=_to_path(cfg["calculate_FND"].get("bg_distance_matrix")),
        output_results=Path(cfg["calculate_FND"].get("output_results", "fnd_results.csv")),
    ) 

    return PipelineConfig(
        inputs=inputs,
        outputs=outputs,
        compute=compute,
        neighbour_divergence=neighbour_div,
        random_neighbour_divergence=random_neighbour_div,
        ortholog_filter=ortholog_filter,
        FND=FND,
        raw=cfg,
    )


def validate_config(c: PipelineConfig) -> None:
    if c.inputs.dist_matrix is None and c.inputs.background_emb is None:
        raise ValueError("Need either inputs.dist_matrix OR inputs.background_emb.")

    if c.inputs.ortholog_emb is None:
        raise ValueError("inputs.ortholog_emb is required.")

    if c.compute.chunk_size <= 0:
        raise ValueError("compute.chunk_size must be > 0.")

    if c.ortholog_filter.enabled:
        path = c.ortholog_filter.homolog_annotation
        if path is None or not path.exists():
            raise ValueError("Invalid ortholog_filter.homolog_annotation: expected a valid file path.")
        if not (0.0 <= c.ortholog_filter.min_species_coverage <= 1.0):
            raise ValueError("ortholog_filter.min_species_coverage must be between 0 and 1.")
        if c.ortholog_filter.min_species < 1:
            raise ValueError("ortholog_filter.min_species must be >= 1.")
        #if c.inputs.homolog_annotation is None:
        #    raise ValueError("ortholog_filter.enabled=True requires inputs.homolog_annotation.")
    if c.random_neighbour_divergence.enabled: 
        if c.random_neighbour_divergence.sample_size is None:
            raise ValueError("random_neighbour_divergence.sample_size must be set when random_neighbour_divergence.enabled=True.")
        if c.random_neighbour_divergence.sample_size < 1:
            raise ValueError("random_neighbour_divergence.sample_size must be >= 1.")
        if c.outputs.out_random is None:
            raise ValueError("random_neighbour_divergence.enabled=True requires outputs.out_random to be set.")
        #if c.random_neighbour_divergence.sample_size > len(c.inputs.background_emb):
        #    raise ValueError("random_neighbour_divergence.sample_size cannot be greater than the number of available genes in background_emb.")
    


def compute_ortholog_divergence(
    engine: ComputeNDistanceDict,
    gene_pos: Iterable[str],
    ortholog_df: pd.DataFrame,
    dist_type: str,
    apply_ortholog_filter: bool,
    return_mode: Literal["aggregate", "per_gene"],
    filter_params: Optional[Dict[str, Any]] = None,
    chunk_size: Optional[int] = None
):
    """Dispatch to within/between ortholog divergence based on cfg.dist_type."""
    
    common_kwargs = dict(
        gene_pos_list=gene_pos,
        ortholog_df=ortholog_df,
        apply_ortholog_filter=apply_ortholog_filter,
        filter_params=filter_params,
    )
    
    if dist_type == "between":
        return engine.between_ortholog_divergence(**common_kwargs,
                                                  chunk_size=chunk_size)

    if dist_type == "within":
        return engine.within_ortholog_divergence(
            **common_kwargs,
            chunk_size = chunk_size,
            group_by_gene=(return_mode == "per_gene"),
        )

    raise ValueError(f"Invalid dist_type: {dist_type}")

##RUN PIPELINES

def run_divergence_pipeline(cfg: PipelineConfig) -> pd.DataFrame:
    """
    Calculate distance matrix 
    """

    # IO setup
    out_dir = Path(cfg.outputs.out_distance)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    
    #LOAD ARGUMENTS
    out_dir = cfg.outputs.out_distance
    dist_type = cfg.neighbour_divergence.dist_type
    segment = cfg.neighbour_divergence.segment
    ortholog_emb_path = cfg.inputs.ortholog_emb
    chunk_size = cfg.compute.chunk_size
    
    targets = load_fam_map(cfg.inputs.fam_map)  # Dict[str, List[str]]
    
    apply_ortholog_filter = cfg.ortholog_filter.enabled
    filter_params = None
    if apply_ortholog_filter:
        homol_ref = read_pickle(cfg.ortholog_filter.homolog_annotation)
        filter_params = dict(
            relationship=cfg.ortholog_filter.ortholog_types,
            min_species_coverage=cfg.ortholog_filter.min_species_coverage,
            min_species=cfg.ortholog_filter.min_species,
            homology_ref = homol_ref
            
        )
    
     # compute state / distance matrix
    bg_embed_df = _load_embeddings(cfg.inputs.background_emb)
    engine = ComputeNDistanceDict(bg_embed_df)
    
    id_to_idx = {sid: i for i, sid in enumerate(engine.ids)}
    
    n_done = 0
    for (fam_id,genes) in targets.items():
        
        present_ids = [g for g in genes if any(g in id for id in id_to_idx)]
        missing = [g for g in genes if not any(g in id for id in id_to_idx)]

        if not present_ids:
            print(f"No gene IDs found for family {fam_id}, skipping.")
            continue
        if missing:
            print(f"Warning: {len(missing)} gene IDs not found for family {fam_id}: {missing}")
        
        ortho_df = _load_ortholog_embeddings(ortho_embeddings_path=ortholog_emb_path, gene_ids=genes)

        grp_ids = _parse_ids(ortho_df.iloc[:,0])
        gene_pos = set(grp_ids.loc[grp_ids["__genes"].isin(genes), "__gpos"])

        # 2) process
        dist = compute_ortholog_divergence(
            engine=engine,
            gene_pos=gene_pos,
            ortholog_df=ortho_df,
            dist_type=dist_type,
            return_mode=cfg.neighbour_divergence.return_mode,
            apply_ortholog_filter=apply_ortholog_filter,
            filter_params= filter_params,
            chunk_size=chunk_size
        )

        # 3) save
        out_path = f"{out_dir}/{fam_id}_{dist_type}_{segment}_diverge"
        save_pickle(dist, out_path)
        print(f"Saved ndist for family {fam_id} to {out_path}")
        n_done += 1

    return pd.DataFrame(rows)


def run_random_divergence_pipeline(cfg:PipelineConfig
) -> pd.DataFrame:
    """
    Calculate distance matrix 
    """
    
    out_dir = cfg.outputs.out_distance
    out_path = cfg.outputs.out_random
    
    dist_type = cfg.random_neighbour_divergence.dist_type
    segment = cfg.random_neighbour_divergence.segment
    sample_size = cfg.random_neighbour_divergence.sample_size
    random_seed = cfg.random_neighbour_divergence.random_seed
    chunk_size = cfg.compute.chunk_size

    #filtering
    apply_ortholog_filter = cfg.ortholog_filter.enabled
    
    #print("Apply ortholog filter?", apply_ortholog_filter)
    print("sample size", sample_size)
    filter_params = None
    if apply_ortholog_filter:
        homol_ref = read_pickle(cfg.ortholog_filter.homolog_annotation)
        filter_params = dict(
            relationship=cfg.ortholog_filter.ortholog_types,
            min_species_coverage=cfg.ortholog_filter.min_species_coverage,
            min_species=cfg.ortholog_filter.min_species,
            homology_ref = homol_ref
        )        
    
    bg_embed_df = _load_embeddings(cfg.inputs.background_emb)
    engine = ComputeNDistanceDict(bg_embed_df)
    
    # Sample genes uniformly
    bg_ids = _parse_ids(engine.ids)
    #genes = bg_ids["__genes"].dropna().unique()
    gpos = bg_ids["__gpos"].dropna().unique()
    if sample_size > len(gpos):
        raise ValueError(f"sample_size={sample_size} > available gene positions={len(gpos)}")

    #rand_genes = (pd.Series(genes).sample(n=sample_size, random_state=random_seed).tolist())
    #rand_gpos = (bg_ids.loc[bg_ids["__genes"].isin(rand_genes), "__gpos"].dropna().unique().tolist())
    
    if random_seed:
        rand_gpos = (pd.Series(gpos).sample(n=sample_size, random_state=random_seed).tolist())
    else:
        rand_gpos = (pd.Series(gpos).sample(n=sample_size).tolist())
    rand_genes = (bg_ids.loc[bg_ids["__gpos"].isin(rand_gpos), "__genes"].dropna().unique().tolist())
    #print('rand_genes', len(rand_genes))
    #load random embeddings
    print('load ortho embeddings...')
    ortho_df= _load_ortholog_embeddings(
                ortho_embeddings_path=cfg.inputs.ortholog_emb,
                gene_ids=rand_genes
            )
    
    # 2) process
    print('compute divergence')
    dist = compute_ortholog_divergence(
                engine=engine,
                gene_pos=rand_gpos,
                ortholog_df=ortho_df,
                dist_type=dist_type,
                return_mode =cfg.random_neighbour_divergence.return_mode,
                apply_ortholog_filter=apply_ortholog_filter,
                filter_params=filter_params,
                chunk_size=chunk_size
            )

    # Save
    if out_path is None:
        out_path = Path(out_dir) / f"random_{sample_size}_{dist_type}_{segment}_diverge"
    save_pickle(dist, str(out_path))
    print(f"Saved to {out_path}.pkl")
    
    return 



def run_FND_pipeline(cfg: "PipelineConfig") -> pd.DataFrame:
    """
    Compute Family Neighbour Divergence (FND), ensuring prerequisite distance outputs exist.

    Inputs
    ------
    cfg.FND.fam_distance_matrix : Optional[Union[str, Path]]
        Family-distance input. May be:
        - a FILE path (single aggregated family distance object), or
        - a DIRECTORY containing per-family distance outputs (e.g. '*diverge.pkl') and a
          sentinel success flag (e.g. '_SUCCESS.json') written when the divergence pipeline completes.

        If None, this function runs `run_divergence_pipeline(cfg)` and uses `cfg.outputs.out_distance`
        as the family-distance directory.

    cfg.FND.bg_distance_matrix : Optional[Union[str, Path]]
        Background/random distance input. Must be a FILE path.
        If None, this function runs `run_random_divergence_pipeline(cfg)` and uses a default output
        file path (recommended to store in config).

    Behavior
    --------
    - If fam_distance_matrix is a directory: waits for a success flag in that directory.
      If it is a file: waits for the file to exist.
    - Always waits for the background distance file to exist.
    - Calls `calc_FND(...)` exactly once.

    Returns
    -------
    pd.DataFrame
        FND results returned by `calc_FND`.
    """
    # Resolve family distance path (dir or file)
    fam_dist_path = cfg.FND.fam_distance_matrix
    if fam_dist_path is None:
        print(f"Distance path for family divergence not provided, generating distances from {cfg.inputs.background_emb}") #and {cfg.inputs.ortholog_emb} using {cfg.inputs.fam_map}")
        run_divergence_pipeline(cfg)
        fam_dist_path = Path(cfg.outputs.out_distance)   # directory expected to contain _SUCCESS.json
    else:
        fam_dist_path = Path(fam_dist_path)

    # Resolve background distance path (file)
    bg_dist_path = cfg.FND.bg_distance_matrix
    if bg_dist_path is None:
        print(f"Distance path for background divergence not provided, running random divergence")
        # Configure + run background/random pipeline
        cfg.random_neighbour_divergence.enabled = True
        cfg.random_neighbour_divergence.sample_size = 100
        cfg.random_neighbour_divergence.dist_type = "within"
        
        # put this in config as cfg.FND.bg_distance_matrix_output or similar
        bg_dist_path = Path(cfg.outputs.out_distance) / "random_bg_distance_matrix.pkl"
        
        cfg.output.out_random =  bg_dist_path

        run_random_divergence_pipeline(cfg)

        
    else:
        bg_dist_path = Path(bg_dist_path)


    # Compute FND
    return calc_FND(
        fam_dist_path,
        bg_dist_path,
        save_path=cfg.FND.output_results,
    )

    

def run_pipeline(cfg_path: str | Path):
    cfg = load_pipeline_config(cfg_path)
    validate_config(cfg)
    
    print(cfg)
    
    if cfg.FND.enabled:
        print('Running Family Neighbour Divergence (FND) calculation ...')
        # For now, just compute FND directly from distance matrix if provided, otherwise compute distance matrix from embeddings and then compute FND, but in the future could add option to compute FND directly from embeddings without computing full distance matrix by using a more efficient algorithm that only computes distances to nearest neighbours instead of all pairwise distances
        run_FND_pipeline(cfg)
        return
    
    #if random divergence, call that pipeline instead
    if cfg.random_neighbour_divergence.enabled:
        print('Running random neighbour divergence...')
        run_random_divergence_pipeline(cfg)
        return
    
    else:
        print('Running family neighbour divergence...')
        # for now just run neighbour divergence, can add CLI args to specify which pipeline to run if needed
        run_divergence_pipeline(cfg)
        return
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run ortholog divergence pipeline.")
    parser.add_argument("config", type=str, help="Path to YAML config file.")
    args = parser.parse_args()

    run_pipeline(args.config)



