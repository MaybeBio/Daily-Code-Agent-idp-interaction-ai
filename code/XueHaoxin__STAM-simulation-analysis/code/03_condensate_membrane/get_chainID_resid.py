import numpy as np

def get_chainID(resindex, tot_chain, single_len, ag):
    """
    Determine the chain ID (A-J) from a residue index.

    Parameters:
        resindex: Zero-based residue index relative to the first protein residue.
        tot_chain: Total number of chains (<=10).
        single_len: Number of residues per chain.
        ag: MDAnalysis AtomGroup used to validate the total residue count.

    Returns:
        str: Chain ID (A, B, C, ..., J).
    """
    total_res = len(ag.residues)
    expect_total = tot_chain * single_len
    if total_res != expect_total:
        raise ValueError(
            f"Residue count mismatch: found {total_res}, expected {expect_total} "
            f"(tot_chain={tot_chain} x single_len={single_len})"
        )
    residues_per_chain = single_len

    chain_index = resindex // residues_per_chain

    chain_mapping = "ABCDEFGHIJ"
    if chain_index >= len(chain_mapping):
        raise ValueError(
            f"Chain index {chain_index} exceeds the A-J mapping range"
        )

    return chain_mapping[chain_index]

def get_local_resid(resindex, tot_chain, single_len):
    """
    Calculate the one-based residue ID within its chain.

    Parameters:
        resindex: Zero-based residue index relative to the first protein residue.
        tot_chain: Retained for compatibility; not used in this calculation.
        single_len: Number of residues per chain.

    Returns:
        int: Residue ID within the chain, starting from 1 for each chain.
    """
    residues_per_chain = single_len
    # Calculate the starting residue index of the corresponding chain
    chain_index = resindex // residues_per_chain  #0-based chain index
    chain_start_resindex = chain_index * residues_per_chain  # Chain start index

    # Calculate the one-based residue ID within the chain
    local_resid = resindex - chain_start_resindex + 1
    return local_resid