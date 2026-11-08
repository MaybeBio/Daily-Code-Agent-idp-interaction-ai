import pandas as pd
import numpy as np
import os
import json
import ast
import xgi
from hypergraphx import Hypergraph

def safe_eval(x):
     try:
          return ast.literal_eval(x) if isinstance(x, str) else x
     except:
          return np.nan


def extract_protein_features(file_path, species=None, species_taxon_id=None):
    # JSON
    if file_path.lower().endswith('.json'):
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Handle nested JSON
        if isinstance(data, dict):
            for key in data.keys():
                if isinstance(data[key], list):
                    data = data[key]
                    break
            else:
                data = [data]

        df = pd.DataFrame(data)

        # Rename
        if 'name' in df.columns:
            df = df.rename(columns={'name': 'protein_name'})

        features = [
            'uniprot_id', 'uniprotkb_id', 'protein_name', 'species_taxon_id',
            'species_name', 'gene_name', 'ensembl_id', 'iupred_score',
            'synthetic_condensate_count'
        ]

        df = df[features]

        # Filter species
        if species is not None and 'species_name' in df.columns:
            df = df[df['species_name'] == species]
        if species_taxon_id is not None and 'species_taxon_id' in df.columns:
            df = df[df['species_taxon_id'] == species_taxon_id]

        # Clean NaNs
        df = df.replace(['n/a', '', '{}'], np.nan)

        # Fix columns
        if 'species_taxon_id' in df.columns:
            df['species_taxon_id'] = df['species_taxon_id'].replace(0, np.nan)

        if 'iupred_score' in df.columns:
            df['iupred_score'] = df['iupred_score'].apply(safe_eval)

        return df

    # CSV
    elif file_path.lower().endswith('.csv'):
        df = pd.read_csv(file_path)

        df = df.rename(columns={
            'Uniprot ID': 'uniprot_id',
            'Gene Name': 'gene_name',
            'Name': 'protein_name',
            'Biomolecular condensate count': 'biomolecular_condensate_count',
            'Synthetic condensate count': 'synthetic_condensate_count',
            'Species': 'species_name'
        })

        if species is not None and 'species_name' in df.columns:
            df = df[df['species_name'] == species]

        return df

    else:
        raise ValueError("Unsupported file format")


def create_protein_dataset(json_path, csv_path, species, species_taxon_id = None, save_path = None):

     if not json_path.lower().endswith('.json'):
          raise ValueError('json_path must be a JSON file')
     if not csv_path.lower().endswith('.csv'):
          raise ValueError('csv_path must be a CSV file')

     protein_json = extract_protein_features(json_path, species, species_taxon_id)
     protein_csv = extract_protein_features(csv_path, species)

     # Merge DataFrames
     df = pd.merge(
          protein_json,
          protein_csv[['uniprot_id']],
          on='uniprot_id',
          how='left'
     )

     if save_path is not None:
          # Extract the directory from the full path
          directory = os.path.dirname(save_path)
          
          # Create directory if it doesn't exist
          if directory:
               os.makedirs(directory, exist_ok=True)
          
          # Save the DataFrame to CSV
          df.to_csv(save_path, index=False)

     return df


def extract_condensate_features(
    file_path,
    confidence_th=None,
    condensate_type=None,
    species=None,
    species_tax_id = None,
    num_proteins=None
):

    if condensate_type is not None and condensate_type not in ['biomolecular', 'synthetic', 'infectious']:
        raise ValueError("condensate_type must be: biomolecular, synthetic, infectious, or None")

    if confidence_th is not None and not (1 <= confidence_th <= 5):
        raise ValueError("confidence_th must be in [1, 5]")

    # JSON
    if file_path.lower().endswith('.json'):
        with open(file_path, 'r') as f:
            data = json.load(f)

        if isinstance(data, dict):
            for key, value in data.items():
                if isinstance(value, list):
                    data = value
                    break
            else:
                data = [data]

        df = pd.DataFrame(data)

        if 'name' in df.columns:
            df = df.rename(columns={'name': 'condensate_name'})

        features = [
            'uid', 'condensate_name', 'species_taxon_id',
            'species_name', 'confidence_score', 'condensate_type'
        ]

        df = df[features]

        # Filters
        if confidence_th is not None:
            df = df[df['confidence_score'] >= confidence_th]

        if condensate_type is not None:
            df = df[df['condensate_type'] == condensate_type]
        
        if species_tax_id is not None:
            df = df[df['species_taxon_id'] == species_tax_id]

        if species is not None:
            df = df[df['species_name'] == species]

        return df

    # CSV
    elif file_path.lower().endswith('.csv'):
        df = pd.read_csv(file_path)

        df = df.rename(columns={
            'Name': 'condensate_name',
            'UID': 'uid',
            'Condensate Type': 'condensate_type',
            'Species Tax Id': 'species_taxon_id',
            'Proteins': 'proteins',
            'Condensatopathy': 'condensatopathy',
            'Confidence Score': 'confidence_score'
        })

        if species_tax_id is not None and 'species_taxon_id' in df.columns:
            df = df[df['species_taxon_id'] == species_tax_id]

        # Filter by number of proteins
        if num_proteins is not None and 'proteins' in df.columns:
            df = df[df['proteins'] >= num_proteins]

        return df

    else:
        raise ValueError("Unsupported file format")


def create_condensate_dataset(
    json_path,
    csv_path,
    confidence_th=None,
    condensate_type=None,
    species=None,
    num_proteins=None,
    species_tax_id = None,
    save_path = None
):

     if not json_path.lower().endswith('.json'):
          raise ValueError('json_path must be a JSON file')
     if not csv_path.lower().endswith('.csv'):
          raise ValueError('csv_path must be a CSV file')

     condensate_json = extract_condensate_features(
          json_path,
          confidence_th,
          condensate_type,
          species,
          species_tax_id
     )

     condensate_csv = extract_condensate_features(
          csv_path,
          confidence_th,
          condensate_type,
          species,
          species_tax_id,
          num_proteins          
     )

     csv_columns = condensate_csv[['uid', 'proteins', 'DNA', 'RNA', 'C-mods', 'condensatopathy']]
     json_columns = condensate_json[['uid', 'species_name', 'confidence_score', 'condensate_type']]
     # Merge DataFrames
     df = pd.merge(json_columns, csv_columns, on='uid', how='left')
     df = df.dropna(subset=['proteins'])

     if save_path is not None:
        # Extract the directory from the full path
        directory = os.path.dirname(save_path)
        
        # Create directory if it doesn't exist
        if directory:
            os.makedirs(directory, exist_ok=True)
        
        # Save the DataFrame to CSV
        df.to_csv(save_path, index=False)

     return df


def create_hypergraph(proteins_df, condensates_df, mapping_df, min_proteins=3, save_path=None):

    proteins_df = proteins_df.copy()
    condensates_df = condensates_df.copy()
    mapping_df = mapping_df.copy()

    # ------------------------------------------------------------------
    # 0. Drop unwanted condensate columns early
    # ------------------------------------------------------------------
    cols_to_drop = ["DNA", "RNA", "C-mods", "condensatopathy"]
    condensates_df = condensates_df.drop(
        columns=[c for c in cols_to_drop if c in condensates_df.columns]
    )

    # ------------------------------------------------------------------
    # 1. Basic cleaning
    # ------------------------------------------------------------------
    mapping_df["uniprotkb_ac"] = mapping_df["uniprotkb_ac"].astype(str).str.strip().str.upper()
    mapping_df["condensate_id"] = mapping_df["condensate_id"].astype(str).str.strip().str.upper()
    proteins_df["uniprot_id"] = proteins_df["uniprot_id"].astype(str).str.strip().str.upper()

    if "uid" in condensates_df.columns:
        condensates_df["uid"] = condensates_df["uid"].astype(str).str.strip().str.upper()

    # Rename proteins to protein_count
    if "proteins" in condensates_df.columns:
        condensates_df = condensates_df.rename(columns={"proteins": "protein_count"})

    if "protein_count" in condensates_df.columns:
        condensates_df = condensates_df.dropna(subset=["protein_count"])
        condensates_df["protein_count"] = condensates_df["protein_count"].astype(int)

    # ------------------------------------------------------------------
    # 2. Merge mapping with condensate metadata
    # ------------------------------------------------------------------
    condensates_mapping = pd.merge(
        mapping_df,
        condensates_df,
        how="left",
        left_on="condensate_id",
        right_on="uid"
    )

    # ------------------------------------------------------------------
    # 3. Merge with protein table
    # ------------------------------------------------------------------
    df = pd.merge(
        condensates_mapping,
        proteins_df,
        how="left",
        left_on="uniprotkb_ac",
        right_on="uniprot_id",
        suffixes=("_cond", "_prot")
    )

    # Clean duplicated species column
    if "species_name_cond" in df.columns:
        df = df.drop(columns=["species_name_cond"])
    if "species_name_prot" in df.columns:
        df = df.rename(columns={"species_name_prot": "species_name"})

    # ------------------------------------------------------------------
    # 4. Fix missing metadata
    # ------------------------------------------------------------------

    # Fill condensate_type for synthetic entries
    if "condensate_type" in df.columns and "condensate_name" in df.columns:
        synthetic_mask = (
            df["condensate_type"].isna() &
            df["condensate_name"].str.contains("Synthetic", case=False, na=False)
        )
        df.loc[synthetic_mask, "condensate_type"] = "synthetic"

    # Ensure uid exists and fill missing values
    if "uid" not in df.columns:
        df["uid"] = df["condensate_id"]
    else:
        df["uid"] = df["uid"].fillna(df["condensate_id"])

    # Keep only matched proteins
    df = df.dropna(subset=["uniprot_id"])

    # ------------------------------------------------------------------
    # 5. Additional cleaning (your new requirements)
    # ------------------------------------------------------------------

    # Fill confidence_score
    if "confidence_score" in df.columns:
        df["confidence_score"] = df["confidence_score"].fillna(0.0)

    # Recompute protein_count from actual data
    if "uniprotkb_ac" in df.columns:
        protein_counts = (
            df.groupby("condensate_id")["uniprotkb_ac"]
            .transform("nunique")
        )
        df["protein_count"] = df["protein_count"].fillna(protein_counts)

    # Drop unwanted columns
    cols_to_drop = ["gene_name", "ensembl_id"]
    df = df.drop(columns=[c for c in cols_to_drop if c in df.columns])

    # Drop rows with missing condensate_type
    if "condensate_type" in df.columns:
        df = df.dropna(subset=["condensate_type"])
    if min_proteins is not None:
        df = df[df['protein_count'] >= min_proteins]

    # ------------------------------------------------------------------
    # 6. Save incidence table
    # ------------------------------------------------------------------
    if save_path is not None:
        directory = os.path.dirname(save_path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        df.to_csv(save_path, index=False)

    # ------------------------------------------------------------------
    # 7. Build hypergraph
    # ------------------------------------------------------------------
    hyperedge_dict = (
        df.groupby("condensate_id")["uniprotkb_ac"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    filtered_edges = {
        edge_id: members
        for edge_id, members in hyperedge_dict.items()
        if len(members) >= min_proteins
    }

    H = xgi.Hypergraph(filtered_edges)

    # ------------------------------------------------------------------
    # 8. Node attributes
    # ------------------------------------------------------------------
    node_cols = ["uniprotkb_ac", "protein_name", "iupred_score"]
    node_cols = [c for c in node_cols if c in df.columns]

    if "uniprotkb_ac" in node_cols:
        node_data = df[node_cols].drop_duplicates(subset=["uniprotkb_ac"]).copy()
        node_attr = node_data.set_index("uniprotkb_ac").to_dict("index")
        node_attr = {n: attrs for n, attrs in node_attr.items() if n in H.nodes}
        H.set_node_attributes(node_attr)

    # ------------------------------------------------------------------
    # 9. Edge attributes
    # ------------------------------------------------------------------
    edge_cols = [
        "condensate_id", "uid", "condensate_name",
        "protein_count", "confidence_score", "condensate_type"
    ]
    edge_cols = [c for c in edge_cols if c in df.columns]

    if "condensate_id" in edge_cols:
        edge_data = df[edge_cols].drop_duplicates(subset=["condensate_id"]).copy()
        edge_data = edge_data[edge_data["condensate_id"].isin(H.edges)]
        edge_attr = edge_data.set_index("condensate_id").to_dict("index")
        H.set_edge_attributes(edge_attr)

    return H


def create_hypergraph_hgx(
    proteins_df,
    condensates_df,
    mapping_df,
    min_proteins=3,
    save_path=None,
):

    proteins_df = proteins_df.copy()
    condensates_df = condensates_df.copy()
    mapping_df = mapping_df.copy()

    # ------------------------------------------------------------
    # 0. Drop unwanted condensate columns
    # ------------------------------------------------------------
    cols_to_drop = ["DNA", "RNA", "C-mods", "condensatopathy"]

    condensates_df = condensates_df.drop(
        columns=[c for c in cols_to_drop if c in condensates_df.columns]
    )

    # ------------------------------------------------------------
    # 1. Basic cleaning
    # ------------------------------------------------------------
    mapping_df["uniprotkb_ac"] = (
        mapping_df["uniprotkb_ac"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    mapping_df["condensate_id"] = (
        mapping_df["condensate_id"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    proteins_df["uniprot_id"] = (
        proteins_df["uniprot_id"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    if "uid" in condensates_df.columns:
        condensates_df["uid"] = (
            condensates_df["uid"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

    # Rename proteins -> protein_count
    if "proteins" in condensates_df.columns:
        condensates_df = condensates_df.rename(
            columns={"proteins": "protein_count"}
        )

    # ------------------------------------------------------------
    # 2. Merge mapping + condensates
    # ------------------------------------------------------------
    condensates_mapping = pd.merge(
        mapping_df,
        condensates_df,
        how="left",
        left_on="condensate_id",
        right_on="uid",
    )

    # ------------------------------------------------------------
    # 3. Merge proteins
    # ------------------------------------------------------------
    df = pd.merge(
        condensates_mapping,
        proteins_df,
        how="left",
        left_on="uniprotkb_ac",
        right_on="uniprot_id",
        suffixes=("_cond", "_prot"),
    )

    # Species cleanup
    if "species_name_cond" in df.columns:
        df = df.drop(columns=["species_name_cond"])

    if "species_name_prot" in df.columns:
        df = df.rename(
            columns={"species_name_prot": "species_name"}
        )

    # ------------------------------------------------------------
    # 4. Metadata fixes
    # ------------------------------------------------------------
    if (
        "condensate_type" in df.columns
        and "condensate_name" in df.columns
    ):
        synthetic_mask = (
            df["condensate_type"].isna()
            & df["condensate_name"].str.contains(
                "Synthetic",
                case=False,
                na=False,
            )
        )

        df.loc[synthetic_mask, "condensate_type"] = "synthetic"

    if "uid" not in df.columns:
        df["uid"] = df["condensate_id"]
    else:
        df["uid"] = df["uid"].fillna(df["condensate_id"])

    df = df.dropna(subset=["uniprot_id"])

    # ------------------------------------------------------------
    # 5. Additional cleaning
    # ------------------------------------------------------------
    if "confidence_score" in df.columns:
        df["confidence_score"] = (
            df["confidence_score"]
            .fillna(0.0)
        )

    protein_counts = (
        df.groupby("condensate_id")["uniprotkb_ac"]
        .transform("nunique")
    )

    df["protein_count"] = protein_counts

    if "condensate_type" in df.columns:
        df = df.dropna(subset=["condensate_type"])

    if min_proteins is not None:
        df = df[df["protein_count"] >= min_proteins]

    # ------------------------------------------------------------
    # 6. Save incidence table
    # ------------------------------------------------------------
    if save_path is not None:
        directory = os.path.dirname(save_path)

        if directory:
            os.makedirs(directory, exist_ok=True)

        df.to_csv(save_path, index=False)

    # ------------------------------------------------------------
    # 7. Build hyperedges
    # ------------------------------------------------------------
    grouped = (
        df.groupby("condensate_id")
    )

    edge_list = []
    edge_metadata = []

    for condensate_id, group in grouped:

        members = tuple(
            sorted(group["uniprotkb_ac"].unique())
        )

        if len(members) < min_proteins:
            continue

        edge_list.append(members)

        meta_cols = [
            "condensate_id",
            "uid",
            "condensate_name",
            "protein_count",
            "confidence_score",
            "condensate_type",
        ]

        meta_cols = [
            c for c in meta_cols if c in group.columns
        ]

        metadata = (
            group[meta_cols]
            .iloc[0]
            .to_dict()
        )

        edge_metadata.append(metadata)

    # ------------------------------------------------------------
    # 8. Node metadata
    # ------------------------------------------------------------
    node_cols = [
        "uniprotkb_ac",
        "protein_name",
        "iupred_score",
    ]

    node_cols = [
        c for c in node_cols if c in df.columns
    ]

    node_metadata = {}

    node_df = (
        df[node_cols]
        .drop_duplicates(subset=["uniprotkb_ac"])
    )

    for _, row in node_df.iterrows():

        node_id = row["uniprotkb_ac"]

        attrs = (
            row.drop(labels=["uniprotkb_ac"])
            .dropna()
            .to_dict()
        )

        node_metadata[node_id] = attrs

    # ------------------------------------------------------------
    # 9. Create HypergraphX hypergraph
    # ------------------------------------------------------------
    H = Hypergraph(
        edge_list=edge_list,
        node_metadata=node_metadata,
        edge_metadata=edge_metadata,
    )

    return H