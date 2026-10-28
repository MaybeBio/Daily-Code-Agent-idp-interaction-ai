from typing import Dict, Iterable, List, Literal, Optional, Union,Any
import numpy as np
import pandas as pd
from scipy import stats
import math
from collections import defaultdict
from itertools import product,combinations
from dataclasses import dataclass
from pathlib import Path
import json
import yaml
import os

from idr_diverge.utils.helpers import transform_data, collect_values_by_key, embedding_h5_to_dataframe, read_pickle, find_files_by_prefix
from idr_diverge.distances.embed_distance import EmbedDistanceMatrix


#Implements the core logic for computing neighbour distances between protein segment embeddings, focuses on comparing orthologous regions across species. Includes family neighbour divergence (FND) calculation.

#This module was written by Denise Le, with some functions edited and refactored with assistance from ChatGPT (OpenAI) for clarity and maintainability. 
#All changes were reviewed and validated by the author.



# ----------------------------
# IO helpers
# ----------------------------


def _parse_ids(ortholog_ids:List[str]|pd.Series):
    """
    Vectorized parse of ID strings into __species and __gpos columns.
    Assumes IDs look like 'SPECIES_..._...'; __gpos is the join of the last |pos| tokens.
    """
    pos= -4
    ids = pd.Series(ortholog_ids).astype(str)  # Series[str]
    parts = ids.str.split('_')                      # Series[list[str]]

    k = abs(pos)
    if k == 0:
        raise ValueError("pos must be non-zero; use negative values to take last k tokens")

    # join the last k tokens (vectorized for list-like via .str.slice)
    gpos = parts.str.slice(-k).str.join('_')
    genes = parts.str.get(-k)
    species_full = parts.str.get(0)                 # first token

    # remove non-letters; if result is empty, fall back to original
    species_clean = species_full.str.replace(r'[^A-Za-z]', '', regex=True)
    species = species_clean.where(species_clean.ne(''), species_full)
    
    df = pd.DataFrame()
    df['__id'] = ids
    df['__species'] = species
    df['__ensembl_id'] = species_full
    df['__genes'] = genes
    df['__gpos'] = gpos
    return df

def chunk_list(lst, chunk_size):
        return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]
    
def _load_embeddings(embed_path: str | Path) -> pd.DataFrame:
    """Load the embeddings (csv or .h5)  as DataFrame from path."""
    path = Path(embed_path)
    if path.suffix == '.csv':
        embed_df = pd.read_csv(path,header=None)
    elif path.suffix == '.h5':
        embed_df = embedding_h5_to_dataframe(path)
    else:
        raise ValueError('Embedding path must be csv or h5')
    return embed_df


def _load_ortholog_embeddings(ortho_embeddings_path: str | Path, gene_ids: Iterable[str]) -> pd.DataFrame:
    """
    Load ortholog embeddings from a file or directory for the given uniprot gene IDs.

    If `embeddings_path` is:
      • a file: loads that file directly
      • a directory: loads only files matching the provided gene IDs (by prefix)

    ParametersFm
    ----------
    embeddings_path : str or Path
        Path to embedding file or directory of embedding files.

    gene_ids : iterable of str
        Gene IDs used to select files when loading from a directory.

    Returns
    -------
    pd.DataFrame
        Concatenated ortholog embeddings.

    Raises
    ------
    FileNotFoundError
        If path does not exist.

    ValueError
        If directory contains no matching files.
    """

    path = Path(ortho_embeddings_path)

    if not path.exists():
        raise FileNotFoundError(f"Ortholog embedding path not found: {path}")

    # Case 1: single embedding file
    if path.is_file():
        return _load_embeddings(path)

    # Case 2: directory of embedding files
    files = find_files_by_prefix(prefixes=gene_ids, input_dir=path)

    if not files:
        raise ValueError(
            f"No embedding files found in directory {path} "
            f"for gene IDs: {list(gene_ids)[:5]}..."
        )

    # Sort ensures deterministic order
    files = sorted(files)
    dfs = [_load_embeddings(f) for f in files]

    return pd.concat(dfs, ignore_index=True)


def load_fam_map(path: Path | str) -> Dict[str, List[str]]:
    path = Path(path)
    suf = path.suffix.lower()

    if suf == ".json":
        obj = json.loads(path.read_text())
        if not isinstance(obj, dict):
            raise ValueError("fam-map JSON must be a dict: {fam_id: [genes...]}")
        return {str(k): [str(x) for x in v] for k, v in obj.items()}

    if suf in [".yml", ".yaml"]:
        if yaml is None:
            raise SystemExit("PyYAML not installed. Run: pip install pyyaml")
        obj = yaml.safe_load(path.read_text())
        if not isinstance(obj, dict):
            raise ValueError("fam-map YAML must be a dict: {fam_id: [genes...]}")
        return {str(k): [str(x) for x in v] for k, v in obj.items()}

    if suf in [".tsv", ".txt", ".csv"]:
        sep = "," if suf == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep)
        required = {"fam_id", "gene"}
        if not required.issubset(df.columns):
            raise ValueError(f"fam-map table must have columns {required}, found: {list(df.columns)}")
        fam_to_genes: Dict[str, List[str]] = {}
        for fam, sub in df.groupby("fam_id"):
            fam_to_genes[str(fam)] = sub["gene"].astype(str).tolist()
        return fam_to_genes

    raise ValueError("Unsupported fam-map format. Use .json/.yml/.tsv/.csv")

def dict_max_depth(obj: Any) -> int:
    """
    Return the maximum depth of nested dictionaries.

    Conventions:
      - Non-dict objects have depth 0
      - An empty dict has depth 1
      - A flat dict (no dict values anywhere inside) has depth 1
    """
    if not isinstance(obj, dict):
        return 0
    if not obj:
        return 1
    child_depths: List[int] = []
    for v in obj.values():
        if isinstance(v, dict):
            child_depths.append(dict_max_depth(v))
        elif isinstance(v, (list, tuple, set)):
            # If dicts are inside iterables, include their depths too
            for item in v:
                if isinstance(item, dict):
                    child_depths.append(dict_max_depth(item))

    return 1 if not child_depths else 1 + max(child_depths)


def _load_distance_dict(dist_path: str) -> Dict[str, Dict[str, List[float]]]:
    """
    Load and normalize a distance dictionary from a pickle file - aggregate distances by species.

    Expected input shapes:
      depth=2: {key -> {...}}  (top-level groups; values are dict-like groups)
      depth=3: {id -> {key -> {...}}}

    Normalization:
      depth=2 => collect_values_by_key(values) - get list of distances aggregated by species
      depth=3 => {id: collect_values_by_key(inner.values())}
    """

    if not os.path.isfile(dist_path):
        raise FileNotFoundError(f"Distance pickle not found: {dist_path}")
    dist_dict = read_pickle(dist_path)
    if not isinstance(dist_dict, dict):
        raise TypeError(f"Expected a dict from {dist_path}, got {type(dist_dict).__name__}")

    depth = dict_max_depth(dist_dict)
    
    if depth == 1:
        return dist_dict

    if depth == 2:
        return collect_values_by_key(list(dist_dict.values()))
    if depth == 3: #distance dictionary for {fam_id: dist_dict}
        return {
            _id: collect_values_by_key(list(inner.values()))
            for _id, inner in dist_dict.items()
            if isinstance(inner, dict)
        }
    raise ValueError(
        f"Unexpected nested dict depth={depth} for {dist_path}. "
        "Expected depth 2 or 3."
    )


# ----------------------------
# Distance computation 
# ----------------------------

def _precompute_background_sorted(dm: np.ndarray, background_idx: np.ndarray) -> list[np.ndarray]:
    """
    Precompute sorted distances from each row to the background embeddings.

    Parameters
    ----------
    distance_matrix : np.ndarray
        Pairwise distance matrix.(cosine)

    background_idx : np.ndarray
        Indices of embeddings that form the background competitor set.

    Returns
    -------
    list[np.ndarray]
        Sorted background distances for each row.
    """
    n = dm.shape[0]
    out = [None] * n
    for i in range(n):
        vals = dm[i, background_idx]
        if np.isnan(vals).any():
            vals = np.where(np.isnan(vals), np.inf, vals)
        out[i] = np.sort(vals, kind="mergesort")
    return out


def _rank_in_row(i:int, j:int, distance:np.ndarray, background_sorted_vals:np.ndarray)->int: 
    """
    Rank position (0-based) of column j in row i under the allowed set:
      allowed = (non-added) ∪ {i if i is added} ∪ {j if j is added}
    Compute by inserting d(i,j) into the pre-sorted base array.
    """
    d = distance[i, j]
    if np.isnan(d):
        return np.inf  
    return int(np.searchsorted(background_sorted_vals[i], d, side='left'))

def neighbour_dist(id1:str, id2:str, distance:np.ndarray, id_to_idx:dict[str,int], background_sorted_vals:np.ndarray) -> float:
    '''
    Neighbour (rank) distance between two embeddings (id1 and id2)
    '''
    i1 = id_to_idx[id1]
    i2 = id_to_idx[id2]
    r12 = _rank_in_row(i1, i2, distance, background_sorted_vals)
    r21 = _rank_in_row(i2, i1, distance, background_sorted_vals)
    return (r12 + r21) / 2.0


def neighbour_dist_groupwise_fast(
    group1_ids: List[str],
    group2_ids: List[str],
    distance_matrix: np.ndarray,
    id_to_idx: dict[str, int],
    background_sorted_vals,  #TODO add type
    return_dict: bool = False,
) -> Union[np.ndarray, dict[tuple[str, str], float]]:
    """
    Compute neighbour distances between all possible pairs between two groups of embeddings.
    Skips pairs where either ID is not present in the ID->index map. 
    Assumes id_to_idx and background_sorted_vals are already precomputed.

    If return_dict=True, returns { (a,b): dist } else returns np.ndarray of distances.
    """

    # Filter to known IDs (dict membership is O(1))
    g1 = [x for x in group1_ids if x in id_to_idx]
    g2 = [x for x in group2_ids if x in id_to_idx]

    if not g1 or not g2:
        return {} if return_dict else np.empty(0, dtype=float)

    if return_dict:
        out: dict[tuple[str, str], float] = {}
        for a, b in product(g1, g2):
            if a == b:
                continue
            out[(a, b)] = neighbour_dist(a, b, distance_matrix, id_to_idx, background_sorted_vals)
        return out

    dists = [
        neighbour_dist(a, b, distance_matrix, id_to_idx, background_sorted_vals)
        for a, b in product(g1, g2)
        if a != b
    ]
    return np.asarray(dists, dtype=float)

# ----------------------------
# Filtering 
# ----------------------------

OrthologRelationship = Literal[
    "ortholog_one2one",
    "ortholog_one2many",
    "ortholog_one2one_one2many",
]

#DEFAULT_HOMOL_REF = read_pickle(CONFIG_DATA['ensembl_homol_ref'])

def filter_ortholog_region_ids(
    ortholog_ids: Iterable[str], #ID format ENSEMBLID_GENE_SEGMENT_start_end
    homology_ref: Dict,
    relationship: OrthologRelationship = "ortholog_one2one",
    min_species_coverage: float = 0.5,
    min_species: int = 10,
    keep_species: str = "HUMAN"
) -> List[str]:
    """
    Filter ortholog region IDs using an ensembl homology reference, per-species coverage, and deduplication.
        
     Parameters
    ----------
    ortholog_ids : Iterable[str]
        Region IDs to filter. ID format:  ENSEMBLID_GENE_SEGMENTTYPE_start_end
    homology_ref : dict
        Ensembl homology reference (ortholog mappings and percent identities).
    relationship : str, optional
        Ortholog type to retain (default: "ortholog_one2one").
    min_species_coverage : float, optional
        Minimum fraction of segments per species.
    min_species : int, optional
        Minimum number of species required.
    keep_species : str, optional
        Species to always retain (default: "HUMAN") Do not modify.
        
    Steps:
      1) Parse region IDs into a DataFrame with columns: __ensembl_id, __species, __gpos, __ids.
      2) Keep rows whose protein is allowed for `relationship`, plus always keep `keep_species`.
      3) Deduplicate within (__species, __gpos): keep highest percent_id; if all missing, pick randomly.
      4) Optionally enforce per-species coverage: each species must cover at least
         ceil(n_unique_gpos * min_species_coverage) unique gene positions.
      5) Require at least `min_species` unique species remaining.
      6) Return filtered original IDs (__ids).
    """

    # ---- 1) Parse ortholog IDs ----
    # Expects columns: __ensembl_id, __species, __gpos, __ids
    df = _parse_ids(ortholog_ids)
    if df.empty:
        return []

    n_unique_gpos = df["__gpos"].nunique()

    # ---- 2) Allowed proteins for orthology relationship ----
    prot_map = homology_ref.get("ensembl_protein_orthologs", {})
    if relationship == "ortholog_one2one_one2many":
        #allowed = set(prot_map.get("ortholog_one2one", [])) | set(prot_map.get("ortholog_one2many", []))
        allowed = set().union(
                                prot_map.get("ortholog_one2one", []),
                                prot_map.get("ortholog_one2many", [])
                            )
    else:
        allowed = set(prot_map.get(relationship, []))

    if not allowed:
        return []

    # Keep allowed proteins, plus always keep keep_species (e.g., HUMAN)
    df = df[df["__ensembl_id"].isin(allowed) | (df["__species"] == keep_species)].copy()
    if df.empty:
        return []

    # ---- 3) Deduplicate within (__species, __gpos) using percent identity ----
    # Check if duplicates exist
    dup_mask = df.duplicated(subset=['__species', '__gpos'], keep=False)

    if dup_mask.any():
        # Attach percent_id for all rows (NaN if missing)
        percent_id_map = homology_ref.get("percent_ids", {})
        df["percent_id"] = df["__ensembl_id"].map(percent_id_map).astype(float)  # missing -> NaN
        #df['percent_id'] = df['__ensembl_id'].apply(lambda x: percent_id_map.get(x, np.nan))
            
        # Randomize row order first so that ties among NaNs (or equal percent_id) are broken randomly
        df = df.sample(frac=1, random_state=42).reset_index(drop=True)

        # Prefer higher percent_id; NaNs go last; keep first row per (__species, __gpos)
        df = (
            df.sort_values(
                by=["__species", "__gpos", "percent_id"],
                ascending=[True, True, False],
                na_position="last",
            )
            .drop_duplicates(subset=["__species", "__gpos"], keep="first")
            .reset_index(drop=True)
        )

    # ---- 4) Enforce species coverage threshold ----
    if min_species_coverage is not None:
        min_required = math.ceil(n_unique_gpos * float(min_species_coverage))
        spp_counts = df.groupby("__species")["__gpos"].nunique() #count unique gpos per species
        #print('spp_counts',spp_counts)
        keep_spp = set(spp_counts[spp_counts >= min_required].index) 
        #print('min_required', min_required)
        #print('keep_spp', keep_spp)
        df = df[df["__species"].isin(keep_spp)].copy()
        if df.empty:
            #print('min species coverage filtered')
            return []

    # ---- 5) Require a minimum number of species ----
    if df["__species"].nunique() < int(min_species):
        #print('not enough species')
        return []

    return df['__id'].tolist()


# ----------------------------
# Core divergence computation (stub)
# ----------------------------
from dataclasses import dataclass

@dataclass
class DistanceState:
    dm: np.ndarray #distance_matrix
    all_ids: list[str]
    id_to_idx: dict[str, int]
    background_sorted_vals: list[np.ndarray]
    background_idx: np.ndarray
    


@dataclass
class ComputeNDistanceDict(): #return distance dictionary
    
    def __init__(self,embed_df,background_ids = None, gpos_tokens =-4, ortholog_filter_params = None):
        embed_df=embed_df.copy()
        self.id_col = embed_df.columns[0]
        self.gpos_tokens = gpos_tokens #what position in the id do you find the gene_position
        
        self.embed_df = embed_df
        #embed_df[id_col] = embed_df[id_col].apply(lambda x: f'HUMAN_{x}').str.upper()       
        self.emb_matrix = EmbedDistanceMatrix(embed_df)
        self.ids =self.emb_matrix.distance_ids

        #background competitor universe for neighbour distance calculation
        if background_ids is None:
            self.background_ids = set(self.ids)   # default = full HUMAN background
        else:
            self.background_ids = set(background_ids)
        
        if ortholog_filter_params is None: #TODO Do I want to keep this
            self.ortholog_filter_params = {'relationship':"ortholog_one2one_one2many", 
                'min_species_coverage':0.3,
                'min_species':5}
    
    def _prepare_extended_dm_state(self, ortholog_df:Optional[pd.DataFrame]=None) -> DistanceState:
        """Update distance matrix and calculate distances relative to the background/base embeddings"""
        emb_matrix = self.emb_matrix
        
        if ortholog_df is None or ortholog_df.empty:
            all_ids = self.ids
            id_to_idx = {id_: k for k, id_ in enumerate(all_ids)}
            
            #Input ids is same as background
            background_idx = np.array(list(id_to_idx.values()),dtype=int)
            
            dm = emb_matrix.distance_matrix
            
            background_sorted_vals = _precompute_background_sorted(dm, background_idx)
            
            return DistanceState(
                dm=dm,
                all_ids=all_ids,
                id_to_idx=id_to_idx,
                background_sorted_vals=background_sorted_vals,
                background_idx = background_idx
            )
        
        
        else:
            ortholog_df = ortholog_df.drop_duplicates(subset=[self.id_col]).copy()
            ortholog_df[self.id_col] = ortholog_df[self.id_col].astype(str).str.strip()

            new_dm, updated_ids = emb_matrix.add_to_distance_matrix(ortholog_df)
            
            new_dm = new_dm.astype(np.float32, copy=False) #Just added 03-24
            
            all_ids = list(updated_ids)
            id_to_idx = {id_: k for k, id_ in enumerate(all_ids)}

            # competitors = background IDs only
            background_idx = np.array(
                [i for i, id_ in enumerate(all_ids) if id_ in self.background_ids],
                dtype=int
            )
            background_sorted_vals = _precompute_background_sorted(new_dm, background_idx)

            return DistanceState(
                dm=new_dm,
                all_ids=all_ids,
                id_to_idx=id_to_idx,
                background_sorted_vals=background_sorted_vals,
                background_idx = background_idx
            )    

    
    
    def between_segment_single_spp(self, gene_pos_list) -> np.ndarray:        
        #distance_matrix = self.emb_matrix.distance_matrix
        ids = self.ids
        select_ids = [id for id in ids if any(g in id for g in gene_pos_list)]
        
        state = self._prepare_extended_dm_state()
        
        n_distances = neighbour_dist_groupwise_fast(
                        group1_ids=select_ids,
                        group2_ids=select_ids,
                        distance_matrix=state.dm,
                        id_to_idx=state.id_to_idx,
                        background_sorted_vals=state.background_sorted_vals,
                        return_dict=False
                    )
        
        return np.array(n_distances)
    
    def between_ortholog_divergence(self,
                         gene_pos_list:List[str],
                         ortholog_df:pd.DataFrame,
                         apply_ortholog_filter: bool = False,
                        filter_params: Optional[dict] = None,
                         chunk_size:int=10) -> Dict:
        id_col =self.id_col

        # Ensure the ID column is string-typed and keep the same order as ids_for_matrix
        ortholog_df[id_col] = ortholog_df[id_col].astype(str).str.strip()
        
        if apply_ortholog_filter: 
            params = filter_params if filter_params is not None else self.ortholog_filter_params

            candidate_ids = ortholog_df[id_col].dropna().unique().tolist()
            filtered_ids = filter_ortholog_region_ids(ortholog_ids = candidate_ids,  **params)

            ortholog_df = ortholog_df[ortholog_df[id_col].isin(filtered_ids)]
            
            if ortholog_df.empty:
                params_str = {k: v for k, v in params.items() if k != 'homology_ref'}
                print(f'No orthologs passed filtering criteria:{params_str}. Gene pos list:{gene_pos_list}')
                return {}       
                                  
        all_ids_combined = list(set(ortholog_df[id_col])|set(self.background_ids)) #TODO change to set(ortholog_df[id_col]) ?
        df = _parse_ids(all_ids_combined)
        
        # Build group map (species, gpos) -> list of IDs
        grp = df.groupby(['__species', '__gpos'])['__id'].apply(list).to_dict()
        all_species = sorted({sp for sp, _ in grp.keys() if sp != 'HUMAN'})
        gene_pos_set = sorted(set(gene_pos_list))
        H = {gp: grp.get(('HUMAN', gp), []) for gp in gene_pos_set}
        S = {gp: [s for sp in all_species for s in grp.get((sp, gp), [])] for gp in gene_pos_set}
        id_to_species = dict(zip(df["__id"], df["__species"]))
        #print(id_to_species)
        
        distances = defaultdict(dict)
        
        #Precompute distances + sort with added ortholog_df                        
        state = self._prepare_extended_dm_state(ortholog_df)
        #orth_ids = ortholog_df[self.id_col].astype(str).str.strip()
        
        combos = list(combinations(gene_pos_set, 2))
        if chunk_size:
            chunks = chunk_list(combos, chunk_size)
            len_combos = len(chunks)
        else: 
            chunks = [combos]
            len_combos = 1
        
        for i, chunk in enumerate(chunks):
            print(f'chunk {i+1}/{len_combos}')
            for (gp1, gp2) in chunk: 
                name = f"{gp1}|{gp2}" #TODO should I change to a tuple?
                #print(name)
                H1, H2 = H[gp1], H[gp2]
                S1, S2 = S[gp1], S[gp2]
                
                if not H1 or not H2:
                    H1 = gp1
                    H2 = gp2
                if not (S1 and S2):#(H1 and H2 and S1 and S2): 
                    continue
                
                # ranks for all pairs         
                ranks1 = neighbour_dist_groupwise_fast(
                        group1_ids=H1,
                        group2_ids=S2,
                        distance_matrix=state.dm, #new_dm
                        id_to_idx=state.id_to_idx,
                        background_sorted_vals=state.background_sorted_vals,
                        return_dict=True
                    )
                
                ranks2 = neighbour_dist_groupwise_fast(
                    group1_ids=H2,
                    group2_ids=S1,
                    distance_matrix=state.dm, #new_dm
                    id_to_idx=state.id_to_idx,
                    background_sorted_vals=state.background_sorted_vals,
                    return_dict=True
                )
                
                #res1= {ortholog: dist for pair, dist in ranks1.items() for ortholog in all_species if ortholog in pair}
                #res2= {ortholog: dist for pair, dist in ranks2.items() for ortholog in all_species if ortholog in pair}
                
                #Get species to distance mapping for each gp pair, taking minimum distance if multiple orthologs per species
                res1_dist = {} 
                for (id1,id2), dist in ranks1.items():
                    ortholog_id = id2 if id_to_species.get(id1) == "HUMAN" else id1
                    spp = id_to_species.get(ortholog_id)

                    if spp and spp != "HUMAN":
                        res1_dist[spp] = min(dist, res1_dist.get(spp, dist))
                
                res2_dist = {}
                for (id1,id2), dist in ranks2.items():
                    ortholog_id = id2 if id_to_species.get(id1) == "HUMAN" else id1
                    spp = id_to_species.get(ortholog_id)

                    if spp and spp != "HUMAN":
                        res2_dist[spp] = min(dist, res2_dist.get(spp, dist))

                agg = defaultdict(list)
                
                for ranks in (res1_dist,res2_dist):
                    for spp,dist in ranks.items():
                        agg[spp].append(dist)
                res=agg
                distances[name] = res
        
        return distances 

    def within_ortholog_divergence(self,
                                    gene_pos_list:List[str],
                                    ortholog_df:pd.DataFrame,
                                        apply_ortholog_filter: bool = False,
                                        filter_params: Optional[dict] = None,
                                    group_by_gene:bool=False, #LEFT OFF HERE 02-03-26
                                    chunk_size:int = 1
                                )-> dict: #TODO: type = distancedict class
            """
        Compute within-ortholog divergence for the given gene positions and ortholog dataframe.

        For each gene position (gp):
        - identify HUMAN IDs (H1) and non-HUMAN ortholog IDs (S1) that share that gp
        - compute neighbor-distance ranks between H1 and S1 using an extended distance-matrix state
        - return either:
            * group_by_gene=True: {gp: {species: [distances...]}}
            * group_by_gene=False: {species: [distances...]} aggregate by species across gene positions

        Notes:
        - ortholog filtering is applied to candidate IDs before distance computation, if enabled.
        - requires: _parse_ids, filter_ortholog_region_ids, neighbour_dist_groupwise_fast,
                    collect_values_by_key, and self._prepare_extended_dm_state().
        """
            id_col=self.id_col
            gene_pos_set = list(set(gene_pos_list))

            # Ensure the ID column is string-typed and keep the same order as ids_for_matrix
            ortholog_df[id_col] = ortholog_df[id_col].astype(str)
            
            if apply_ortholog_filter: 
                params = filter_params if filter_params is not None else self.ortholog_filter_params

                candidate_ids = ortholog_df[id_col].dropna().unique().tolist()
                #print(candidate_ids[:10])
                filtered_ids = filter_ortholog_region_ids(ortholog_ids = candidate_ids,  **params)
                #print(filtered_ids[:10])
                ortholog_df = ortholog_df[ortholog_df[id_col].isin(filtered_ids)]
                
                if ortholog_df.empty:
                    param_dict = {k: v for k, v in params.items() if k != 'homology_ref'}
                    print(f'No orthologs passed filtering criteria:{param_dict}. Gene pos list:{gene_pos_set}')
                    return {} 
            
            df = _parse_ids(ortholog_df[id_col])
            
            # Build group map (species, gpos) -> list of IDs
            grp = df.groupby(['__species', '__gpos'])['__id'].apply(list).to_dict()
            all_species = sorted({sp for sp, _ in grp.keys() if sp != 'HUMAN'})
            #gene_pos_set = set(gene_pos_list)
            H = {gp: grp.get(('HUMAN', gp), []) for gp in gene_pos_set}
            #print('H:',H)
            S = {gp: [s for sp in all_species for s in grp.get((sp, gp), [])] for gp in gene_pos_set}
            id_to_species = dict(zip(df["__id"], df["__species"]))
            
            species_dist_list= []
            species_dist_by_gpos = defaultdict(dict)
            
            #gene_pos_set = set(df['__gpos'].unique()) & set(gene_pos_set)
            gene_pos_set = set(gene_pos_set)
            
            #print('gene pos list:',gene_pos_set)
            #print('ortholog df gpos:', set(df['__gpos'].unique()))
            
            if not gene_pos_set:
                return ValueError(f"Ortholog_df ids do not overlap with gene_pos_set: {gene_pos_set}")
            
            if chunk_size:
                chunks = chunk_list(list(gene_pos_set), chunk_size)
                len_combos = len(chunks)
            else: 
                chunks = list(list(gene_pos_set))
                len_combos = 1
            
            for i, gp_chunk in enumerate(chunks):
                state=None 
                print(f'chunk {i+1}/{len_combos}')
                #print(gp_chunk)
                chunk_ids = df.loc[df["__gpos"].isin(gp_chunk), "__id"].unique()
                #print('Chunk ids:', gp_chunk, chunk_ids[:10])
                ortholog_chunk_df = ortholog_df[ortholog_df[id_col].isin(chunk_ids)]
                
                if ortholog_chunk_df.empty:
                    print('empty')
                    continue
                
                #Precompute distances + sort with added ortholog_df   
                print('prepare extended dm')
                                    
                state = self._prepare_extended_dm_state(ortholog_chunk_df)
                
                for gp in gp_chunk:            
                    H1 = H[gp] #human id for this gene pos
                    S1 = S[gp] #species orthologs for this gene pos
                    
                    if not H1:
                        H1 = gp
                        #print(f'skipping {gp} because no human orthologs')
                    elif not S1:
                        print(f'skipping {gp} because no species orthologs')
                        continue
            
                    # ranks for all pairs
                    ranks = neighbour_dist_groupwise_fast(
                            group1_ids=H1,
                            group2_ids=S1,
                            distance_matrix=state.dm, #new_dm
                            id_to_idx=state.id_to_idx,
                            background_sorted_vals=state.background_sorted_vals,
                            return_dict=True
                        )
                    
                    #Get species to distance mapping for each gp pair, taking minimum distance if multiple orthologs per species
                    species_dist = {} #human to ortholog species distance
                    for (id1,id2), dist in ranks.items():
                        ortholog_id = id2 if id_to_species.get(id1) == "HUMAN" else id1
                        spp = id_to_species.get(ortholog_id)

                        if spp and spp != "HUMAN":
                            #handle duplicates by taking minimum distance for each species (should not have duplicates after filtering, but just in case)
                            species_dist[spp] = min(dist, species_dist.get(spp, dist))
                        
                    if group_by_gene:
                        species_dist_by_gpos[gp] = species_dist 
                    else:
                        species_dist_list.append(species_dist)
                
            if group_by_gene:
                return species_dist_by_gpos
            else: #aggregate by species
                species_dist_by_species  = collect_values_by_key(species_dist_list)
                return species_dist_by_species

        


# ----------------------------
# Family neighbour divergence computation
# ----------------------------

def differences(fam_dist: Dict[str, np.ndarray], #TODO is is np array or float?
                bg_dist: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """
    Compute element-wise differences between matching entries in two dictionaries.
    For each shared key (ortholog pair), returns:
        fam_dist[key] -  bg_dist[key]
    Keys with missing or empty values are skipped.

    Parameters
    ----------
    fam_dist : dict
        Dictionary of primary distance/divergence values.
    bg_dist : dict
        Dictionary of comparison (background) values.

    Returns
    -------
    Dictionary of differences for shared keys (ortholog pairs).
    """
    
    diff: Dict[str, np.ndarray] = {}

    common_keys = fam_dist.keys() & bg_dist.keys()

    for ortho in common_keys:
        fam_val = fam_dist.get(ortho, ())
        bg_val  = bg_dist.get(ortho, ())

        if fam_val.size == 0 or bg_val.size == 0:
            continue

        diff[ortho] = fam_val - bg_val

    return diff

def _load_fam_distance_dict(path: str | Path, endswith = ".pkl") -> dict:

    path = Path(path)

    if path.is_dir():
        files = sorted(path.glob(f"*{endswith}"))

        if not files:
            raise FileNotFoundError(f"No '*{endswith}' files found in {path}")

        return {
            p.name.split("_", 1)[0]: _load_distance_dict(p)
            for p in files
        }

    else:
        return {path.name.split("_", 1)[0]:_load_distance_dict(path)}

def calc_FND(fam_dist_dict_path,
            unrelated_dist_dict_path,
            data_transform:Literal['log','log_mean','geo_mean','mean', 'median','log_median','geo_median'] = 'log_mean',
            save_path=None) -> pd.DataFrame: 
    
    '''
    Compute Family Neighbour Divergence (FND) for each protein family (in fam_dist_dict).

    FND is defined as the mean difference between transformed neighbour distances 
    of between paralogous segments (between cross-species paralogs) and a background distribution 
    (within orthologs of unrelated segments), computed per ortholog and averaged across each family.

    Parameters
    ----------
    fam_dist_dict_path : str or Path
        Path to family distance dictionary.
    unrelated_dist_dict_path : str or Path
        Path to background distance dictionary.
    data_transform : str, optional
        Transformation applied to distances before comparison (default: 'log_mean').
    save_path : str or Path, optional
        If provided, appends results to CSV.

    Returns
    -------
    pd.DataFrame
        One row per family with:
        - Group_id
        - per-ortholog divergence values (dict)
        - FND (mean divergence)
        - one-sample t-test statistic and p-value
    '''
    
    fam_dist_dict_path = Path(fam_dist_dict_path)
    unrelated_dist_dict_path = Path(unrelated_dist_dict_path)
    
    if not fam_dist_dict_path.exists():
        raise FileNotFoundError(f"Family distance dictionary not found at {fam_dist_dict_path}")
    if not unrelated_dist_dict_path.exists():
        raise FileNotFoundError(f"Background distance dictionary not found at {unrelated_dist_dict_path}")
    
    unrelated_dist_dict = _load_distance_dict(unrelated_dist_dict_path)
    unrelated_dist_dict= transform_data(unrelated_dist_dict,type_=data_transform)

    fam_dist_dict = _load_fam_distance_dict(fam_dist_dict_path, endswith=".pkl")
    
    rows = []
    for fam_id, fam_dists in fam_dist_dict.items():
        if fam_dists is None:
            continue
        #print(fam_id)
        fam_summ_data = transform_data(fam_dists, type_=data_transform)
        diff = differences(fam_summ_data, unrelated_dist_dict)
    
        row = {
            "Group_id": fam_id,
            f"{data_transform} Divergence_per_ortholog": diff,
        }

        if diff:
            vals = np.asarray(list(diff.values()), dtype=float)
            row[f"FND"] = float(vals.mean())
            t_stat, p_val = stats.ttest_1samp(vals, popmean=0)
            row["1samp_ttest_stat"] = float(t_stat)
            row["1samp_ttest_pvalue"] = float(p_val)

        rows.append(row)
        
        if save_path:
            # Convert to DataFrame and save to CSV
            family_divergence_df = pd.DataFrame([row])
            family_divergence_df.to_csv(save_path, index=False, mode='a', header=not os.path.exists(save_path))
            family_divergence_df = None

    return pd.DataFrame(rows)