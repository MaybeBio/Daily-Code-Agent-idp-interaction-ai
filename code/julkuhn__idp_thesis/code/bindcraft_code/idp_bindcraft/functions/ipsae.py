"""ipSAE and related interface-confidence scores.

Reimplementation of the metrics from Dunbrack, "Estimating the accuracy of
protein-protein interactions with AlphaFold and other prediction programs"
(bioRxiv 2025.02.10.637595), validated numerically against the reference
ipsae.py (DunbrackLab/IPSAE, v4).

Why this differs from the project's existing i_pTM_iface: that metric picks the
interface by a 4 A all-atom distance cutoff and uses a single d0 from the total
interface size. ipSAE picks pairs by PAE alone (no distance), and its primary
variant normalises d0 per residue by that residue's own number of confident
partners. The point of the paper is that this removes the arbitrary distance
cutoff and stops a large flexible chain from diluting the score -- the same
problem the i_*_iface metrics were introduced to work around here.

Chain order convention throughout: AF2 binder protocol emits target first, then
binder, and `chains` must follow the same residue order as `pae_matrix`.
"""

import math

import numpy as np


def calc_d0(L, min_value=1.0):
    """TM-score d0 normalisation for an alignment of length L."""
    L = float(L)
    if L > 27:
        d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    else:
        d0 = 1.0
    return max(min_value, d0)


def calc_d0_array(L, min_value=1.0):
    """Vectorised d0. Mirrors the reference, which clamps L at 26 before the cube root."""
    L = np.maximum(26.0, np.asarray(L, dtype=float))
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def ptm_func(x, d0):
    """pTM kernel: 1/(1+(PAE/d0)^2)."""
    return 1.0 / (1.0 + (np.asarray(x, dtype=float) / d0) ** 2.0)


def _asym_scores(pae_matrix, chains, chain1, chain2, pae_cutoff, d0_floor=1.0):
    """All ipSAE variants for the directed pair chain1 -> chain2.

    Aligning on residues of chain1, scoring residues of chain2. Each variant
    differs only in which d0 normalises the pTM kernel. `d0_floor` is 1.0 for
    protein pairs and 2.0 when either chain is nucleic acid.
    """
    chains = np.asarray(chains)
    numres = len(chains)

    in1 = chains == chain1
    in2 = chains == chain2

    # Pair selection is by PAE only -- no distance cutoff anywhere in ipSAE.
    valid_pairs_matrix = np.outer(in1, in2) & (pae_matrix < pae_cutoff)

    # d0chn: normalise by the full size of both chains (this is ipTM-like).
    n0chn = int(in1.sum() + in2.sum())
    d0chn = calc_d0(n0chn, d0_floor)
    ptm_d0chn = ptm_func(pae_matrix, d0chn)

    # d0dom: normalise by the residues that actually take part in confident pairs.
    res1 = np.unique(np.where(valid_pairs_matrix.any(axis=1) & in1)[0])
    res2 = np.unique(np.where(valid_pairs_matrix.any(axis=0) & in2)[0])
    n0dom = int(len(res1) + len(res2))
    d0dom = calc_d0(n0dom, d0_floor)
    ptm_d0dom = ptm_func(pae_matrix, d0dom)

    # d0res: per-residue normalisation by that residue's own partner count.
    n0res_byres = valid_pairs_matrix.sum(axis=1)
    d0res_byres = calc_d0_array(n0res_byres, d0_floor)

    iptm_d0chn_byres = np.zeros(numres)
    ipsae_d0chn_byres = np.zeros(numres)
    ipsae_d0dom_byres = np.zeros(numres)
    ipsae_d0res_byres = np.zeros(numres)

    for i in range(numres):
        if not in1[i]:
            continue
        valid = valid_pairs_matrix[i]
        # ipTM analogue: every residue of chain2, regardless of PAE.
        if in2.any():
            iptm_d0chn_byres[i] = ptm_d0chn[i, in2].mean()
        if valid.any():
            ipsae_d0chn_byres[i] = ptm_d0chn[i, valid].mean()
            ipsae_d0dom_byres[i] = ptm_d0dom[i, valid].mean()
            ipsae_d0res_byres[i] = ptm_func(pae_matrix[i], d0res_byres[i])[valid].mean()

    # The reported asymmetric score is the best-scoring aligned residue, not a mean.
    def _best(byres):
        if not in1.any():
            return 0.0, None
        idx = int(np.argmax(np.where(in1, byres, -np.inf)))
        return float(byres[idx]), idx

    iptm_d0chn, iptm_res = _best(iptm_d0chn_byres)
    ipsae_d0chn, d0chn_res = _best(ipsae_d0chn_byres)
    ipsae_d0dom, d0dom_res = _best(ipsae_d0dom_byres)
    ipsae_d0res, d0res_res = _best(ipsae_d0res_byres)

    return {
        "ipsae": ipsae_d0res,
        "ipsae_d0dom": ipsae_d0dom,
        "ipsae_d0chn": ipsae_d0chn,
        "iptm_d0chn": iptm_d0chn,
        "n0chn": n0chn,
        "n0dom": n0dom,
        "n0res": int(n0res_byres[d0res_res]) if d0res_res is not None else 0,
        "d0chn": d0chn,
        "d0dom": d0dom,
        "d0res": float(d0res_byres[d0res_res]) if d0res_res is not None else 0.0,
        "ipsae_residue": d0res_res,
        "n_valid_pairs": int(valid_pairs_matrix.sum()),
        "byres": {
            "ipsae": ipsae_d0res_byres,
            "ipsae_d0dom": ipsae_d0dom_byres,
            "ipsae_d0chn": ipsae_d0chn_byres,
            "iptm_d0chn": iptm_d0chn_byres,
        },
    }


def compute_ipsae(pae_matrix, chains, chain1="A", chain2="B", pae_cutoff=10.0,
                  d0_floor=1.0):
    """ipSAE for one chain pair, in both directions plus the reported maximum.

    pae_matrix : (N, N) PAE in Angstroms, same residue order as `chains`.
    chains     : length-N sequence of chain labels.

    Returns a flat dict; `ipsae` is the headline score (d0res variant, max over
    the two directions), matching the reference script's `ipsae_max`.
    """
    pae_matrix = np.asarray(pae_matrix, dtype=float)
    chains = np.asarray(chains)

    if pae_matrix.shape[0] != len(chains):
        raise ValueError(
            f"PAE matrix ({pae_matrix.shape[0]}) and chain labels ({len(chains)}) disagree"
        )

    fwd = _asym_scores(pae_matrix, chains, chain1, chain2, pae_cutoff, d0_floor)
    rev = _asym_scores(pae_matrix, chains, chain2, chain1, pae_cutoff, d0_floor)

    out = {}
    for key in ("ipsae", "ipsae_d0dom", "ipsae_d0chn", "iptm_d0chn"):
        out[key] = max(fwd[key], rev[key])
        out[f"{key}_asym_{chain1}{chain2}"] = fwd[key]
        out[f"{key}_asym_{chain2}{chain1}"] = rev[key]

    # Report the geometry of whichever direction won the headline score.
    best = fwd if fwd["ipsae"] >= rev["ipsae"] else rev
    for key in ("n0chn", "n0dom", "n0res", "d0chn", "d0dom", "d0res", "n_valid_pairs"):
        out[key] = best[key]

    out["_byres"] = {chain1: fwd["byres"], chain2: rev["byres"]}
    out["_asym"] = {(chain1, chain2): fwd, (chain2, chain1): rev}
    return out


# Residues the reference treats as nucleic acid; a chain containing any of them
# gets a d0 floor of 2.0 instead of 1.0.
NUC_RESIDUES = frozenset(
    ["DA", "DC", "DG", "DT", "DU", "DI", "A", "C", "G", "U", "I", "N"]
)


def chain_types(chains, residue_types):
    """'nucleic_acid' for any chain containing a nucleotide residue, else 'protein'."""
    chains = np.asarray(chains)
    residue_types = np.asarray(residue_types)
    out = {}
    for chain in np.unique(chains):
        res = residue_types[chains == chain]
        out[str(chain)] = ("nucleic_acid"
                           if any(r in NUC_RESIDUES for r in res) else "protein")
    return out


def pair_d0_floor(chain_types_map, chain1, chain2):
    """d0 floor for a chain pair: 2.0 if either chain is nucleic acid, else 1.0."""
    if (chain_types_map.get(chain1) == "nucleic_acid"
            or chain_types_map.get(chain2) == "nucleic_acid"):
        return 2.0
    return 1.0


def compute_all_pairs(pae_matrix, chains, pae_cutoff=10.0, residue_types=None):
    """ipSAE for every unordered chain pair in the structure.

    Returns {(chain1, chain2): result} keyed by the sorted pair, where each
    result is what compute_ipsae returns (both directions plus the max).
    Nucleic-acid chains get the reference's d0 floor of 2.0 when
    `residue_types` is supplied.
    """
    chains = np.asarray(chains)
    uniq = [str(c) for c in np.unique(chains)]
    ctypes = chain_types(chains, residue_types) if residue_types is not None else {}

    out = {}
    for i, a in enumerate(uniq):
        for b in uniq[i + 1:]:
            out[(a, b)] = compute_ipsae(
                pae_matrix, chains, a, b, pae_cutoff=pae_cutoff,
                d0_floor=pair_d0_floor(ctypes, a, b),
            )
    return out


def compute_pdockq(pae_matrix, chains, coords, plddt, chain1="A", chain2="B",
                   dist_cutoff=8.0):
    """pDockQ and pDockQ2, the metrics ipSAE is benchmarked against.

    coords : (N, 3) CB coordinates (CA for glycine), same order as `chains`.
    plddt  : (N,) per-residue pLDDT on the 0-100 scale.

    Unlike ipSAE these do use a distance cutoff (8 A on CB).
    """
    pae_matrix = np.asarray(pae_matrix, dtype=float)
    chains = np.asarray(chains)
    coords = np.asarray(coords, dtype=float)
    plddt = np.asarray(plddt, dtype=float)

    in1 = chains == chain1
    in2 = chains == chain2

    dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    contact = np.outer(in1, in2) & (dist <= dist_cutoff)

    npairs = int(contact.sum())
    if npairs == 0:
        # keep the key set identical to the populated return, so callers that read
        # pdockq2_asym do not have to special-case non-contacting chain pairs
        return {"pdockq": 0.0, "pdockq2": 0.0, "pdockq2_asym": [0.0, 0.0],
                "n_contact_pairs": 0}

    iface = np.zeros(len(chains), dtype=bool)
    iface |= contact.any(axis=1)
    iface |= contact.any(axis=0)
    mean_plddt = float(plddt[iface].mean())

    x = mean_plddt * math.log10(npairs)
    pdockq = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018

    # pDockQ2 replaces the pair count with the mean pTM kernel at a fixed d0 of 10.
    # PAE is asymmetric, so this differs by direction; the reference reports the max.
    pdockq2_asym = []
    for block in (pae_matrix[contact], pae_matrix.T[contact]):
        mean_ptm = float(ptm_func(block, 10.0).mean())
        x2 = mean_plddt * mean_ptm
        pdockq2_asym.append(1.31 / (1 + math.exp(-0.075 * (x2 - 84.733))) + 0.005)

    return {
        "pdockq": pdockq,
        "pdockq2": max(pdockq2_asym),
        "pdockq2_asym": pdockq2_asym,
        "n_contact_pairs": npairs,
    }


def compute_lis(pae_matrix, chains, chain1="A", chain2="B", cutoff=12.0):
    """LIS (local interaction score): mean of (12-PAE)/12 over sub-cutoff interchain pairs.

    The reference reports the average of the two directions, not the maximum.
    """
    pae_matrix = np.asarray(pae_matrix, dtype=float)
    chains = np.asarray(chains)

    scores = []
    for a, b in ((chain1, chain2), (chain2, chain1)):
        block = pae_matrix[np.ix_(chains == a, chains == b)]
        valid = block[block < cutoff]
        scores.append(float(((cutoff - valid) / cutoff).mean()) if valid.size else 0.0)

    return {"lis": float(np.mean(scores)), "lis_asym": scores}
