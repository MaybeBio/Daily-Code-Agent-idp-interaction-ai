from collections import defaultdict
from typing import Iterable, Dict, Any, List
from pathlib import Path
import pandas as pd
import h5py
import re
from typing import List, Dict, Any, Iterable, Union
import pickle



def load_config(config_path):
    '''Read config yaml as dictionary.'''
    import yaml
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def resolve_config_paths(config, config_path, path_keys=None):
    """
    Recursively resolve relative path values in a config dictionary.

    Parameters
    ----------
    config : dict
        Parsed YAML config.
    config_path : str or Path
        Path to the YAML file.
    path_keys : set[str] | None
        Optional set of keys whose values should be treated as paths.
        If None, path-like strings are inferred automatically.

    Returns
    -------
    dict
        New config dictionary with resolved paths.
    """
    config_dir = Path(config_path).resolve().parent

    def is_pathlike(key, value):
        if not isinstance(value, str):
            return False

        # Explicit key list takes priority
        if path_keys is not None:
            return key in path_keys

        # Otherwise infer from value
        p = Path(value)
        return (
            "/" in value
            or "\\" in value
            or value.startswith(".")
            or p.suffix != ""
        )

    def recurse(obj, parent_key=None):
        if isinstance(obj, dict):
            return {
                k: recurse(v, k)
                for k, v in obj.items()
            }

        if isinstance(obj, list):
            return [recurse(v, parent_key) for v in obj]

        if is_pathlike(parent_key, obj):
            p = Path(obj)
            if not p.is_absolute():
                return config_dir / p
            return p

        return obj

    return recurse(config)


##      Read and write files    ##


#Read hdf5 files with 'ids' and 'embeddings' datasets and return a DataFrame with one row per sequence and embedding columns. Optionally change the ID based on the filename.
def embedding_h5_to_dataframe(file_path: str|Path, change_id=False): #TODO remove change_id option -> should preprocess embedding files instead of preprocessing each time
    """
    Loads an HDF5 file with 'ids' and 'embeddings' datasets and returns a DataFrame.

    Parameters:
    - file_path (str): Path to the .h5 file.
    - change_id (bool): Whether to modify the ID based on the filename.

    Returns:
    - pd.DataFrame: DataFrame with one row per sequence and embedding columns.
    """
    
    file_path = str(file_path)  # Ensure it's a string for h5py
    
    with h5py.File(file_path, "r") as f:
        ids = f["ids"][()]
        embeddings = f["embeddings"][()]

    # Decode byte strings if necessary
    if isinstance(ids[0], bytes):
        ids = [id.decode('utf-8') for id in ids]

    # Optionally change IDs based on the filename
    if change_id:
        filename = file_path.split('/')[-1]
        gene = filename.split('_')[0]
        pos = '_'.join(filename.split('_')[-4:-1])
        ids = [f"{id.split('_')[0]}_{gene}_{pos}" for id in ids]

    # Create DataFrame with IDs and embeddings
    df = pd.DataFrame({
        "id": ids,
        "embedding": list(embeddings)  # ensure it's a list of arrays
    })
    
    # paired_data = list(zip(ids, embeddings)) 
    # df = pd.DataFrame(paired_data, columns=["id", "embedding"])

    # Expand the embedding column into separate numeric columns
    id_col = df["id"]
    df = pd.DataFrame(df["embedding"].tolist())
    
    # Combine ID column as column 0
    df_expanded = pd.concat([id_col, df], axis=1)
    df_expanded.rename(columns={"id": 0}, inplace=True)
    df_expanded.columns = range(df_expanded.shape[1])

    return df_expanded



def save_pickle(data,outpath):
    with open(f"{outpath}.pkl", "wb") as f:
        pickle.dump(data, f)
    

def read_pickle(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    
    return data



def extract_IDs_from_end(
    ids: List[str],
    split_char: str = "_",
    select: Union[str, Iterable[str]]  = ("spp", "gene", "gene_pos"),
    spp_letters_only: bool = True,
) -> Dict[str, Any]:
    """
    Parse IDs by anchoring to the trailing pattern SEGMENT_START_END, then
    infer gene (last token before suffix) and species (everything before gene).
    Example IDs:
      - HOMO_SAPIENS_TP53_idr_10_50
      - TP53_idr_10_50
      - ENSCHYP00000007622_Q07002_IDR_1_93
    """
    if isinstance(select, str):
        select = [select]
    # suffix: segment (no underscores) + start + end
    suf = re.compile(rf"(?P<segment>[^{re.escape(split_char)}]+)"
                     rf"{re.escape(split_char)}(?P<start>-?\d+)"
                     rf"{re.escape(split_char)}(?P<end>-?\d+)$")

    supported = {"spp", "gene", "segment", "start", "end", "gene_pos", "full"}
    if set(select) - supported:
        raise ValueError("Unsupported field in `select`.")

    parsed, invalid = [], []

    for raw in map(str.strip, ids):
        m = suf.search(raw)
        if not m:
            invalid.append(raw); continue

        segment = m.group("segment")
        try:
            start, end = int(m.group("start")), int(m.group("end"))
        except ValueError:
            invalid.append(raw); continue

        prefix = raw[:m.start()].rstrip(split_char)
        if not prefix:  # no gene present
            invalid.append(raw); continue

        tokens = prefix.split(split_char)
        gene = tokens[-1]
        spp_full = split_char.join(tokens[:-1]) or None

        # BUG in your version: you filtered *prefix* (species+gene).
        # Filter only species, optionally letters-only.
        if spp_full and spp_letters_only:
            spp = re.sub(r"[^A-Za-z]", "", spp_full)
            if not spp:  # if letters-only collapses to empty, keep original
                spp = spp_full
        else:
            spp = spp_full

        gene_pos = f"{gene}{split_char}{segment}{split_char}{start}{split_char}{end}"

        parsed.append({
            "full": raw, "spp": spp, "gene": gene,
            "segment": segment, "start": start, "end": end,
            "gene_pos": gene_pos,
        })

    # uniques for requested fields
    unique: Dict[str, List[Any]] = {}
    for field in select:
        vals = {row[field] for row in parsed}
        unique[field] = sorted(vals) if field in {"start", "end"} \
                        else sorted(vals, key=lambda x: (x is None, str(x)))

    return unique


from pathlib import Path
def find_files_by_prefix(prefixes:List[str], input_dir:str, endswith: str=None) -> List[str]:
    """
    Select files from a directory (including subdirectories) whose filenames start with any of the given prefixes.

    Args:
        prefixes (list of str): List of filename prefixes to match (e.g., ['P12345', 'Q99999']).
        input_dir (str): Path to the root directory to search in (searched recursively).
        endswith (str): End of file name (ex. '.fa')

    Returns:
        list of str: Paths to files whose base filenames start with any of the given prefixes.
    """
    input_dir = Path(input_dir)

    if endswith:
        candidates = input_dir.rglob(f"*{endswith}")
    else:
        candidates = input_dir.rglob("*")

    return [str(path)
        for path in candidates
        if path.is_file() and any(path.name.startswith(p) for p in prefixes)
    ]


def find_files_by_identifier(root_path: str, select_ids: List[str], by: str = 'gene') -> List[str]:
#select_paths(root_path: str, select_ids: list, by: str = 'gene'):
        import glob
        """
        Select file paths from a directory based on matching criteria.

        Parameters
        ----------
        root_path : str
            The directory containing the files.
        select_ids : list
            List of IDs to match against file names.
        by : str
            Criterion for matching. Must be one of:
            - 'gene'      : match by gene name
            - 'gene_pos'  : match by gene name + position
            - 'spp'       : match by species prefix

        Returns
        -------
        list
            List of file paths matching the given criteria.
        """
        valid_by = ('gene', 'gene_pos', 'spp')
        if by not in valid_by:
            raise ValueError(f"'by' must be one of {valid_by}, got '{by}'.")

        # Extract the relevant IDs based on the selection mode
        ids = extract_IDs_from_end(ids=select_ids, select=by)[by]
        paths = []
        if by == 'gene_pos':
            for gene_pos in ids:
                gene = gene_pos.split('_')[0]
                pos = '_'.join(gene_pos.split('_')[1:])
                paths.extend(glob.glob(f"{root_path}/{gene}*{pos}*"))
        elif by == 'gene':
            for gene in ids:
                paths.extend(glob.glob(f"{root_path}/*{gene}*"))
        elif by == 'spp':
            for spp in ids:
                paths.extend(glob.glob(f"{root_path}/*{spp}*"))

        return paths



def transform_data(data_dict:Dict, type_:str = 'log_mean'):
    import numpy as np
    
    """
    Transforms a dictionary of lists based on the specified type.
    
    Supported types:
    - 'log': apply log10(x + 1) to each value
    - 'log_mean': apply log10(x + 1), then take the mean
    - 'geo_mean': apply log10(x + 1), then take 10 ** mean
    - 'mean': take the mean of each list
    """
    
    if type_ == 'log':
        return {k: np.log10([x + 1 for x in v]) for k, v in data_dict.items()}
    
    elif type_ == 'log_mean':
        return {k: np.mean(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'geo_mean':
        return {k: 10 ** np.mean(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'mean':
        return {k: np.mean(v) for k, v in data_dict.items()}
    
    elif type_ == 'median':
        return {k: np.median(v) for k, v in data_dict.items()}
    
    elif type_ == 'log_median':
        return {k: np.median(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'geo_median':
        return {k: 10 ** np.median(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    
    else:
        raise ValueError(f"Unsupported type: {type_}. Must be one of {'log','log_mean','geo_mean','mean', 'median','log_median','geo_median'}")
    
    


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary with protein IDs as keys and sequences as values.

    Parameters:
    file_path (str): The path to the FASTA file.

    Returns:
    dict: A dictionary where keys are protein IDs and values are sequences.
    """
    fasta_dict = {}
    with open(file_path, 'r') as file:
        protein_id = None
        sequence_lines = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if protein_id is not None:
                    fasta_dict[protein_id] = ''.join(sequence_lines)
                protein_id = line[1:].split()[0]  # Extract the protein ID
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Don't forget to add the last protein to the dictionary
        if protein_id is not None:
            fasta_dict[protein_id] = ''.join(sequence_lines)
    
    return fasta_dict


#Functions for preprocessing:

def write_h5_embed_from_df(embed_df, out_path):
    import h5py
    import numpy as np
    ids = embed_df.iloc[:, 0].astype(str).to_numpy()
    embeddings = embed_df.iloc[:, 1:].to_numpy(dtype=np.float32)
    
    with h5py.File(out_path, "w") as f:
        # variable-length UTF-8 strings
        str_dtype = h5py.string_dtype(encoding="utf-8")

        f.create_dataset("ids", data=ids, dtype=str_dtype)
        f.create_dataset("embeddings", data=embeddings)
    
    return out_path


def seq_identity(aln1,aln2):

    # Calculate sequence identity
    sequence_identity = sum(a == b for a, b in zip(aln1, aln2)) / len(aln1)

    return sequence_identity

# def combine_dicts(dict_list):
#     '''Combine dictionary of dictionaries by common keys'''
#     combined = defaultdict(list)
#     for d in dict_list:
#         for k, v in d.items():
#             if isinstance(v,list):
#                 combined[k].extend(v)  # append the list
#             else:
#                 combined[k].extend([v])
#     return dict(combined)


def collect_values_by_key(dicts: Iterable[Dict[Any, Any]]) -> Dict[Any, List[Any]]:
    """
    Aggregate values from multiple dictionaries by key.

    For each key appearing in one or more dictionaries, all associated values
    are collected into a list. If a value is itself a list, its elements are
    extended into the result list.

    This function preserves all values and does not overwrite data.

    Example:
        >>> dicts = [{"A": 1, "B": 2},{"A": 4, "C": 5}]
        >>> collect_values_by_key(dicts)
        {'A': [1, 4], 'B': [2], 'C': [5]}

    Args:
        dicts: A list of dictionaries with potentially overlapping keys.

    Returns:
        A dictionary mapping each key to a list of all values encountered.
    """
    aggregated = defaultdict(list)

    for d in dicts:
        for key, value in d.items():
            if isinstance(value, list):
                aggregated[key].extend(value)
            else:
                aggregated[key].append(value)

    return dict(aggregated)


def merge_dict(dicts: Iterable[Dict[Any, Any]]) -> Dict[Any, Any]:
    """
    Merge multiple dictionaries into one, overwriting values on key collisions.

    When the same key appears in multiple dictionaries, the value from the
    last dictionary in the iterable is retained
      Example:
        >>> dicts = [
        ...     {"A": 1, "B": 2},
        ...     {"A": 3},
        ...     {"C": 4}
        ... ]
        >>> merge_dicts_last_wins(dicts)
        {'A': 3, 'B': 2, 'C': 4}
    Args:
        dict_ls: A list of dictionaries to merge.
    Returns:
        A single dictionary containing all key-value pairs from the input dictionaries.
    """
    merged_dict ={}
    for d in dicts:
        merged_dict.update(d)
    return merged_dict
