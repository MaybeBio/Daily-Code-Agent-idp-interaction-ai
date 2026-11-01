#!/usr/bin/env python3

import os
import pandas as pd


def fusion_csv(directory, gene):
    """
    Automatically merged two CSV files in a folder.

    Parameters:
        directory (str): Path to the folder containing the CSV files
        gene (str): Name of the gene (e.g., "TBL1X")
    Returns:
        Merged DataFrame
    """
    
    #  List of CSV files
    csv_files = [f for f in os.listdir(directory) if f.endswith(".csv")]
    
    path1 = os.path.join(directory, csv_files[0])
    path2 = os.path.join(directory, csv_files[1])

    df1 = pd.read_csv(path1, index_col=0, sep="\t", quotechar="\"", encoding="utf-8")
    df2 = pd.read_csv(path2, index_col=0, sep="\t", quotechar="\"", encoding="utf-8")

    # Concatenation
    df = pd.concat([df1, df2], axis=0)

    # Save
    df.to_csv(f"{directory}BioGRID_merged_{gene}.csv", sep="\t")
    print(f"Merged CSV file saved in: {directory}BioGRID_merged_{gene}.csv")
    
    return df


def merged_csv_ods(fusion, ods, gene):
    """
    Merged the BioGRID CSV with the ODS file containing additional data.
    
    Parameters:
        fusion (DataFrame): The merged DataFrame from the fusion_csv function
        ods (str): Path to the ODS file containing additional data
        gene (str): Name of the gene (e.g., "TBL1X")
    Returns:
        Merged DataFrame with ODS data
    """
    
    # List of files
    df = fusion.reset_index()
    ods_file = pd.read_excel(ods, engine="odf")
    ods_file = ods_file.set_index("Name")
    
    # Normalize the interactor names in the BioGRID DataFrame
    df["Official Symbol Interactor A"] = df["Official Symbol Interactor A"].str.strip().str.upper()
    df["Official Symbol Interactor B"] = df["Official Symbol Interactor B"].str.strip().str.upper()
    
    for col in ods_file.columns:
        # Mapping with Interactor A
        df[col] = df["Official Symbol Interactor A"].map(ods_file[col])
        
        # Add Interactor B if NA
        df[col] = df[col].fillna(df["Official Symbol Interactor B"].map(ods_file[col]))

    # Save
    df.to_csv(f"final_merged_{gene}.csv", sep="\t", index=False)
    print(f"Merged CSV file with ODS data saved in: final_merged_{gene}.csv")


def nbr_interactor(df_A, df_B, clef, gene_A, gene_B):
    """
    Function to create a new column in df, merged and count the unique occurrences in a key column
    
    Parameters:
        df_A (DataFrame): The first DataFrame
        df_B (DataFrame): The second DataFrame
        clef (str): Name of the key column for which to count unique values
        gene_A (str): Name of the first gene (e.g., "TBL1X")
        gene_B (str): Name of the second gene (e.g., "TBL1XR1")
    Returns:        
        n_unique (int): Number of unique values in the key column after merging the two CSV files
    """
    
    if df_A and df_B and isinstance(df_A, str) and isinstance(df_B, str):
        df_A = pd.read_csv(df_A, sep="\t", quotechar="\"", encoding="utf-8")
        df_B = pd.read_csv(df_B, sep="\t", quotechar="\"", encoding="utf-8")
    
    # Create the new column by combining Interactor A and Interactor B and merge the two DataFrames
    df_A[clef] = df_A.apply(lambda row: f"{row["Official Symbol Interactor A"]}-{row["Official Symbol Interactor B"]}" if row["Official Symbol Interactor A"] < row["Official Symbol Interactor B"] else f"{row["Official Symbol Interactor B"]}-{row["Official Symbol Interactor A"]}", axis=1)
    df_B[clef] = df_B.apply(lambda row: f"{row["Official Symbol Interactor A"]}-{row["Official Symbol Interactor B"]}" if row["Official Symbol Interactor A"] < row["Official Symbol Interactor B"] else f"{row["Official Symbol Interactor B"]}-{row["Official Symbol Interactor A"]}", axis=1)
    df = pd.concat([df_A, df_B], ignore_index=True)
    
    # Remove the redundant entries from the key column
    df[clef] = df[clef].str.upper()
    if gene_B == "":
        df[clef] = df[clef].str.replace(f"{gene_A}-{gene_A}", f"{gene_A}", regex=False)
        df[clef] = df[clef].str.replace(f"-{gene_A}", "", regex=False)
        df[clef] = df[clef].str.replace(f"{gene_A}-", "", regex=False)
    else:
        df[clef] = df[clef].str.replace(f"{gene_A}-{gene_B}" or f"{gene_B}-{gene_A}", f"{gene_A}{gene_B}", regex=False)
        df[clef] = df[clef].str.replace(f"{gene_A}-{gene_A}", f"{gene_A}", regex=False)
        df[clef] = df[clef].str.replace(f"{gene_B}-{gene_B}", f"{gene_B}", regex=False)
        df[clef] = df[clef].str.replace(f"-{gene_B}", "", regex=False)
        df[clef] = df[clef].str.replace(f"-{gene_A}", "", regex=False)
        df[clef] = df[clef].str.replace(f"{gene_B}-", "", regex=False)
        df[clef] = df[clef].str.replace(f"{gene_A}-", "", regex=False)

    # Count the unique values in the key column
    n_unique = df[clef].nunique()
    
    # Save the merged DataFrame to a new CSV file
    if gene_B == "":
        df.to_csv(f"merged_interactors_{gene_A}.csv", sep="\t", index=False)
    else:
        df.to_csv(f"merged_interactors_{gene_A}_{gene_B}.csv", sep="\t", index=False)
    
    return n_unique  


# Main
gene_A = "TBL1X"
gene_B = "TBL1XR1"
#ods = "partners_TBL1&R.ods"
clef = "Official Symbol Interactor A-B"

fusion_A = fusion_csv(f"{gene_A}/", gene_A)
fusion_B = fusion_csv(f"{gene_B}/", gene_B)

#merged_csv_ods(fusion_A, ods, gene_A)
#merged_csv_ods(fusion_B, ods, gene_B)

interactor = nbr_interactor(f"BioGRID_merged_{gene_A}.csv", f"BioGRID_merged_{gene_B}.csv", clef, gene_A, gene_B)
print(f"Number of unique interactors with '{gene_A}' and '{gene_B}' in human and mouse = {interactor}")

gene = "NFIL3"
interactor = nbr_interactor(f"{gene}/BIOGRID-GENE-110855-5.0.257.tab3.csv", f"{gene}/BIOGRID-GENE-201749-5.0.257.tab3.csv", clef, gene, "")
print(f"Number of unique interactors with {gene} in human and mouse = {interactor}")
