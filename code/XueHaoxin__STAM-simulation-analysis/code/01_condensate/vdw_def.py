import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS, NSResults
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_chainID(resindex, tot_chain, u):
    """
    Return the chain ID (A-J) for a global residue index.
    Parameters:
        resindex: zero-based residue index
        tot_chain: Total number of chains (up to 10).
        u: MDAnalysis Universe
 	Returns:
        str: Chain ID (A、B、C...J)
    """
    
    # Calculate residues per chain and check divisibility
    total_res = len(u.residues)
    if total_res % tot_chain != 0:
        raise ValueError(
            f"Total residue count {total_res} is not divisible by chain count {tot_chain}. "
            f"Check tot_chain or the residue-to-chain assignment."
        )
    
    residues_per_chain = total_res // tot_chain
    
    chain_index = resindex // residues_per_chain
    if chain_index >= tot_chain:
        chain_index = tot_chain - 1  
    
    # Map the chain index to a letter
    chain_mapping = "ABCDEFGHIJ"
    if chain_index >= len(chain_mapping):
        raise ValueError(f"Chain index {chain_index} exceeds the A-J mapping range")
    
    return chain_mapping[chain_index]

def get_local_resid(resindex, tot_chain, u):
    """
    Convert a residue index to a one-based within-chain residue index.
    Parameters:
        resindex: zero-based residue index 
        tot_chain: Total number of chains.
        u: MDAnalysis Universe
    Returns:
        int: 1-based residue index within its chain.
    """
    total_res = len(u.residues)
    if total_res % tot_chain != 0:
        raise ValueError(f"Total residue count {total_res} is not divisible by chain count {tot_chain}")
    residues_per_chain = total_res // tot_chain
    
    chain_index = resindex // residues_per_chain  
    chain_start_resindex = chain_index * residues_per_chain  
    
    local_resid = resindex - chain_start_resindex + 1
    return local_resid


def aggregate_by_local_resid(contact_matrix, tot_chain,len_chain, u):
    """
    Sum a residue contact matrix by within-chain residue indices.
    
    Parameters:
        contact_matrix (np.ndarray): Input square matrix of contact counts or weights.
        tot_chain (int): Total number of chains.
        u: MDAnalysis Universe
    
    Returns:
        agg_matrix (np.ndarray): Aggregated matrix with shape (len_chain, len_chain).
    """
    
    agg_matrix = np.zeros((len_chain, len_chain), dtype=float)
    
    nres = contact_matrix.shape[0]
    
    for i in range(nres):
        for j in range(i, nres):

            contact_value = contact_matrix[i, j]
            
            if contact_value == 0:
                continue
                
            local_resid_i = get_local_resid(i, tot_chain, u)
            local_resid_j = get_local_resid(j, tot_chain, u)
            
            idx_i = local_resid_i - 1
            idx_j = local_resid_j - 1
            
            if local_resid_i == local_resid_j:
                agg_matrix[idx_i, idx_j] += contact_value
            else:
                agg_matrix[idx_i, idx_j] += contact_value
                agg_matrix[idx_j, idx_i] += contact_value

    return agg_matrix
