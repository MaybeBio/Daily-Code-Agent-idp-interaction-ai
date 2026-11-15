#!/usr/bin/env python
"""Command-line ipSAE scoring for AlphaFold2 / AlphaFold3 / Boltz output files.

Wraps the validated metrics in ipsae.py with the file handling and the three
output files of the reference implementation (DunbrackLab/IPSAE), so results are
directly diffable against it.

Usage:
    python ipsae_cli.py <pae.json>  <model.pdb>  <pae_cutoff> <dist_cutoff>   # AF2
    python ipsae_cli.py <pae.json>  <model.cif>  <pae_cutoff> <dist_cutoff>   # AF3
    python ipsae_cli.py <pae.npz>   <model.cif>  <pae_cutoff> <dist_cutoff>   # Boltz
    python ipsae_cli.py <pae.npz>   <model.pdb>  <pae_cutoff> <dist_cutoff>   # Boltz

Writes, beside the structure file:
    <stem>_<pae>_<dist>.txt         one summary row per directed chain pair, plus a max row
    <stem>_<pae>_<dist>_byres.txt   per-residue scores, for locating the interacting region
    <stem>_<pae>_<dist>.pml         PyMOL script with one colouring alias per chain pair

Only the AF2 path is numerically validated against the reference here (that is
the format the test data uses); the AF3 and Boltz readers follow the reference's
logic but should be spot-checked against it before being relied on.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipsae import (  # noqa: E402
    calc_d0, chain_types, compute_all_pairs, compute_lis, compute_pdockq,
    pair_d0_floor, ptm_func,
)

CHAIN_COLOR = {
    "A": "magenta", "B": "marine", "C": "lime", "D": "orange",
    "E": "yellow", "F": "cyan", "G": "lightorange", "H": "pink",
    "I": "deepteal", "J": "forest", "K": "lightblue", "L": "slate",
    "M": "violet", "N": "arsenic", "O": "iodine", "P": "silver",
    "Q": "red", "R": "sulfur", "S": "purple", "T": "olive",
    "U": "palegreen", "V": "green", "W": "blue", "X": "palecyan",
    "Y": "limon", "Z": "chocolate",
}


# --------------------------------------------------------------------- parsing

def _residue_rows(atoms):
    """Collapse atom records into per-residue rows, preserving file order."""
    order, rows = [], {}
    for a in atoms:
        key = (a["chain"], a["resnum"])
        if key not in rows:
            rows[key] = {"chain": a["chain"], "resnum": a["resnum"],
                         "res": a["resname"], "ca": None, "cb": None,
                         "plddt": 0.0, "ca_atom": None, "cb_atom": None}
            order.append(key)
        r = rows[key]
        name = a["atom"]
        if name == "CA" or (r["ca"] is None and name == "C3'"):
            r["ca"] = a["xyz"]
            r["plddt"] = a["bfac"]
            r["ca_atom"] = a["index"]
        # the reference uses CB for distances, falling back to CA for glycine,
        # and C3' for nucleic acids
        if name == "CB" or "C3" in name or (r["res"] == "GLY" and name == "CA"):
            r["cb"] = a["xyz"]
            r["cb_atom"] = a["index"]
    out = [rows[k] for k in order]
    for r in out:
        if r["ca"] is None:
            raise ValueError(f"residue {r['chain']}{r['resnum']} has no CA/C3' atom")
        if r["cb"] is None:
            r["cb"] = r["ca"]
            r["cb_atom"] = r["ca_atom"]
    return out


def read_pdb(path):
    atoms = []
    idx = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            try:
                bfac = float(line[60:66])
            except ValueError:
                bfac = 0.0
            atoms.append({
                "atom": line[12:16].strip(), "resname": line[17:20].strip(),
                "chain": line[21].strip() or "A", "resnum": int(line[22:26]),
                "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                "bfac": bfac, "index": idx,
            })
            idx += 1
    if not atoms:
        raise ValueError(f"no ATOM records in {path}")
    return _residue_rows(atoms)


def read_cif(path):
    """Minimal mmCIF _atom_site reader (AF3 / Boltz output)."""
    fields, atoms, in_loop = {}, [], False
    idx = 0
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                fields[s.split(".", 1)[1]] = len(fields)
                in_loop = True
                continue
            if in_loop:
                if s.startswith("#") or s.startswith("loop_") or not s:
                    if atoms:
                        break
                    in_loop = False
                    continue
                if not (s.startswith("ATOM") or s.startswith("HETATM")):
                    continue
                parts = s.split()
                if len(parts) < len(fields):
                    continue

                def get(*names, default=None):
                    for nm in names:
                        if nm in fields:
                            return parts[fields[nm]]
                    return default

                chain = get("auth_asym_id", "label_asym_id", default="A")
                resnum = get("auth_seq_id", "label_seq_id", default="0")
                try:
                    resnum = int(resnum)
                except ValueError:
                    continue
                try:
                    bfac = float(get("B_iso_or_equiv", default="0") or 0)
                except ValueError:
                    bfac = 0.0
                atoms.append({
                    "atom": get("label_atom_id", "auth_atom_id", default="").strip('"'),
                    "resname": get("label_comp_id", "auth_comp_id", default=""),
                    "chain": chain, "resnum": resnum,
                    "xyz": [float(get("Cartn_x")), float(get("Cartn_y")),
                            float(get("Cartn_z"))],
                    "bfac": bfac, "index": idx,
                })
                idx += 1
    if not atoms:
        raise ValueError(f"no _atom_site records in {path}")
    return _residue_rows(atoms)


def read_pae(path, residues, nres):
    """PAE matrix (and pLDDT where the PAE file carries it) for `residues`."""
    if path.endswith(".npz"):
        data = np.load(path)
        key = "pae" if "pae" in data else list(data.keys())[0]
        pae = np.array(data[key], dtype=float)
        plddt = None
        # Boltz writes pLDDT to a sibling file with 'pae' swapped for 'plddt'
        sib = path.replace("pae", "plddt")
        if sib != path and os.path.exists(sib):
            raw = np.load(sib)
            vals = np.array(raw[list(raw.keys())[0]], dtype=float)
            if vals.max() <= 1.0:
                vals = 100.0 * vals
            plddt = vals
        return pae, plddt

    with open(path) as fh:
        data = json.load(fh)
    if isinstance(data, list):           # some AF2 pipelines wrap it in a list
        data = data[0]

    if "pae" in data:
        pae = np.array(data["pae"], dtype=float)
    elif "predicted_aligned_error" in data:
        pae = np.array(data["predicted_aligned_error"], dtype=float)
    else:
        raise ValueError(f"no 'pae' or 'predicted_aligned_error' key in {path}")

    plddt = None
    if "plddt" in data:                  # AF2: already per-residue
        plddt = np.array(data["plddt"], dtype=float)
    elif "atom_plddts" in data:          # AF3: per-atom, pick the CB (CA for Gly) atom
        atom_plddts = np.array(data["atom_plddts"], dtype=float)
        try:
            plddt = np.array([atom_plddts[r["cb_atom"]] for r in residues])
        except IndexError:
            plddt = None
    return pae, plddt


def align_pae(pae, nres, residues):
    """Reduce a per-token AF3/Boltz PAE to one row/column per residue."""
    if pae.shape[0] == nres:
        return pae
    if pae.shape[0] > nres:
        # tokens beyond the standard residues (ligands, multi-token residues) are
        # dropped by selecting the CA token of each residue
        idx = [r["ca_atom"] for r in residues]
        if max(idx) < pae.shape[0] and len(set(idx)) == nres:
            return pae[np.ix_(idx, idx)]
    raise ValueError(
        f"PAE matrix is {pae.shape[0]}x{pae.shape[0]} but the structure has {nres} "
        f"residues, and they could not be reconciled"
    )


# ------------------------------------------------------------------ formatting

def contiguous_ranges(numbers):
    """'3-7+11+20-22' style residue selection for PyMOL."""
    nums = sorted(set(int(n) for n in numbers))
    if not nums:
        return ""
    out, start, end = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == end + 1:
            end = n
        else:
            out.append(f"{start}" if start == end else f"{start}-{end}")
            start = end = n
    out.append(f"{start}" if start == end else f"{start}-{end}")
    return "+".join(out)


SUMMARY_HEADER = ("\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  "
                  "ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  "
                  "n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   "
                  "dist2  Model\n")

BYRES_HEADER = ("i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT  "
                "    n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  "
                "ipSAE_d0chn ipSAE_d0dom    ipSAE \n")


def summary_row(c1, c2, kind, s, iptm_af, pdockq, pdockq2, lis,
                nres1, nres2, dist1, dist2, pae_s, dist_s, stem):
    return (f'{c1}    {c2}     {pae_s:3}  {dist_s:3}  {kind:5} '
            f'{s["ipsae"]:8.6f}    {s["ipsae_d0chn"]:8.6f}    {s["ipsae_d0dom"]:8.6f}    '
            f'{iptm_af:5.3f}    {s["iptm_d0chn"]:8.6f}    '
            f'{pdockq:8.4f}   {pdockq2:8.4f}   {lis:8.4f}   '
            f'{int(s["n0res"]):5d}  {int(s["n0chn"]):5d}  {int(s["n0dom"]):5d}  '
            f'{s["d0res"]:6.2f}  {s["d0chn"]:6.2f}  {s["d0dom"]:6.2f}  '
            f'{nres1:5d}   {nres2:5d}   {dist1:5d}   {dist2:5d}   {stem}\n')


# ------------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="ipSAE scoring for AF2/AF3/Boltz outputs")
    ap.add_argument("pae_file", help="PAE json (AF2/AF3) or npz (Boltz)")
    ap.add_argument("structure", help="model .pdb or .cif")
    ap.add_argument("pae_cutoff", type=float, nargs="?", default=10.0)
    ap.add_argument("dist_cutoff", type=float, nargs="?", default=15.0)
    ap.add_argument("--iptm", type=float, default=-1.0,
                    help="the platform's own ipTM, reported verbatim in the ipTM_af column")
    ap.add_argument("--outdir", default=None,
                    help="where to write the three output files (default: beside the structure)")
    args = ap.parse_args(argv)

    if args.structure.endswith(".pdb"):
        residues = read_pdb(args.structure)
        stem_base = args.structure[: -len(".pdb")]
    elif args.structure.endswith(".cif"):
        residues = read_cif(args.structure)
        stem_base = args.structure[: -len(".cif")]
    else:
        ap.error("structure must be a .pdb or .cif file")

    nres = len(residues)
    pae, plddt_from_pae = read_pae(args.pae_file, residues, nres)
    pae = align_pae(pae, nres, residues)

    chains = np.array([r["chain"] for r in residues])
    restypes = np.array([r["res"] for r in residues])
    coords = np.array([r["cb"] for r in residues], dtype=float)
    plddt = (np.asarray(plddt_from_pae, dtype=float)
             if plddt_from_pae is not None and len(plddt_from_pae) == nres
             else np.array([r["plddt"] for r in residues], dtype=float))

    pae_s = f"{int(args.pae_cutoff):02d}"
    dist_s = f"{int(args.dist_cutoff):02d}"
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.structure))
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.join(outdir, os.path.basename(stem_base))
    path_stem = f"{stem}_{pae_s}_{dist_s}"
    pdb_stem = os.path.basename(stem_base)

    results = compute_all_pairs(pae, chains, args.pae_cutoff, residue_types=restypes)
    ctypes = chain_types(chains, restypes)

    with open(f"{path_stem}.txt", "w") as OUT, \
            open(f"{path_stem}_byres.txt", "w") as OUT2, \
            open(f"{path_stem}.pml", "w") as PML:
        OUT.write(SUMMARY_HEADER)
        PML.write("# " + SUMMARY_HEADER.lstrip("\n"))
        OUT2.write(BYRES_HEADER)

        dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

        for (ca, cb_) in sorted(results):
            res = results[(ca, cb_)]
            floor = pair_d0_floor(ctypes, ca, cb_)
            pdq = compute_pdockq(pae, chains, coords, plddt, ca, cb_,
                                 dist_cutoff=8.0)
            lis = compute_lis(pae, chains, ca, cb_)

            # per-CHAIN involvement counts, maxed over both directions -- this is
            # what the reference puts in the max row's nres/dist columns
            chain_res = {ca: [], cb_: []}
            chain_dist = {ca: [], cb_: []}

            for (c1, c2) in ((ca, cb_), (cb_, ca)):
                s = res["_asym"][(c1, c2)]
                in1, in2 = chains == c1, chains == c2
                vp = (np.outer(in1, in2) & (pae < args.pae_cutoff))
                nres1 = int((vp.any(axis=1) & in1).sum())
                nres2 = int((vp.any(axis=0) & in2).sum())
                # the dist columns additionally require a CB contact
                vd = vp & (dmat < args.dist_cutoff)
                dist1 = int((vd.any(axis=1) & in1).sum())
                dist2 = int((vd.any(axis=0) & in2).sum())
                chain_res[c1].append(nres1)
                chain_res[c2].append(nres2)
                chain_dist[c1].append(dist1)
                chain_dist[c2].append(dist2)

                i_lis = lis["lis_asym"][0 if (c1, c2) == (ca, cb_) else 1]
                i_pdq2 = pdq["pdockq2_asym"][0 if (c1, c2) == (ca, cb_) else 1]
                row = summary_row(c1, c2, "asym", s, args.iptm, pdq["pdockq"],
                                  i_pdq2, i_lis, nres1, nres2, dist1, dist2,
                                  pae_s, dist_s, pdb_stem)
                OUT.write(row)
                PML.write("# " + row)

                # PyMOL alias: colour each side's confident interface residues
                sel1 = contiguous_ranges(
                    [residues[i]["resnum"] for i in np.where(vp.any(axis=1) & in1)[0]])
                sel2 = contiguous_ranges(
                    [residues[i]["resnum"] for i in np.where(vp.any(axis=0) & in2)[0]])
                # Deliberate deviation from the reference: when a chain pair has no
                # confident residues the reference emits "resi None", which is not
                # valid PyMOL and makes the alias error out. "resi 0" selects
                # nothing, which is what was meant.
                PML.write(
                    f'alias color_{c1}_{c2}, color gray80, all; '
                    f'color {CHAIN_COLOR.get(c1, "magenta")}, chain  {c1} and resi {sel1 or "0"}; '
                    f'color {CHAIN_COLOR.get(c2, "marine")}, chain  {c2} and resi {sel2 or "0"}\n\n')

            # The max row, once per unordered pair. The reference maximises each
            # column independently, so n0res/d0res come from whichever direction
            # won ipSAE and n0dom/d0dom from whichever won ipSAE_d0dom -- these
            # need not be the same direction.
            fwd, rev = res["_asym"][(ca, cb_)], res["_asym"][(cb_, ca)]
            win_res = fwd if fwd["ipsae"] >= rev["ipsae"] else rev
            win_dom = fwd if fwd["ipsae_d0dom"] >= rev["ipsae_d0dom"] else rev
            maxrow = summary_row(
                ca, cb_, "max",
                {"ipsae": res["ipsae"], "ipsae_d0chn": res["ipsae_d0chn"],
                 "ipsae_d0dom": res["ipsae_d0dom"], "iptm_d0chn": res["iptm_d0chn"],
                 "n0res": win_res["n0res"], "d0res": win_res["d0res"],
                 "n0dom": win_dom["n0dom"], "d0dom": win_dom["d0dom"],
                 "n0chn": fwd["n0chn"], "d0chn": fwd["d0chn"]},
                args.iptm, pdq["pdockq"], pdq["pdockq2"], lis["lis"],
                max(chain_res[ca]), max(chain_res[cb_]),
                max(chain_dist[ca]), max(chain_dist[cb_]),
                pae_s, dist_s, pdb_stem)
            OUT.write(maxrow)
            PML.write("# " + maxrow)
            OUT.write("\n")

        # The per-residue file is grouped by ordered chain pair (aligned chain
        # outer, scored chain inner), which is a different traversal from the
        # summary above, so it gets its own pass.
        uniq = sorted({str(c) for c in chains})
        for c1 in uniq:
            for c2 in uniq:
                if c1 == c2:
                    continue
                res = results[tuple(sorted((c1, c2)))]
                s = res["_asym"][(c1, c2)]
                byres = s["byres"]
                floor = pair_d0_floor(ctypes, c1, c2)
                vp = (np.outer(chains == c1, chains == c2) & (pae < args.pae_cutoff))
                for i in range(nres):
                    if chains[i] != c1:
                        continue
                    n0res_i = int(vp[i].sum())
                    d0res_i = calc_d0(max(n0res_i, 26), floor) if n0res_i else floor
                    OUT2.write(
                        f'{i+1:<4d}    {c1:4}      {c2:4}      '
                        f'{residues[i]["resnum"]:4d}           {residues[i]["res"]:3}        '
                        f'{plddt[i]:8.2f}         '
                        f'{int(s["n0chn"]):5d}  {int(s["n0dom"]):5d}  {n0res_i:5d}  '
                        f'{s["d0chn"]:8.3f}  {s["d0dom"]:8.3f}  {d0res_i:8.3f}   '
                        f'{byres["iptm_d0chn"][i]:8.4f}    {byres["ipsae_d0chn"][i]:8.4f}    '
                        f'{byres["ipsae_d0dom"][i]:8.4f}    {byres["ipsae"][i]:8.4f}\n')

    print(f"wrote {path_stem}.txt")
    print(f"wrote {path_stem}_byres.txt")
    print(f"wrote {path_stem}.pml")
    for (ca, cb_), res in sorted(results.items()):
        print(f"  {ca}-{cb_}: ipSAE={res['ipsae']:.6f} "
              f"(A->B {res[f'ipsae_asym_{ca}{cb_}']:.6f}, "
              f"B->A {res[f'ipsae_asym_{cb_}{ca}']:.6f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
