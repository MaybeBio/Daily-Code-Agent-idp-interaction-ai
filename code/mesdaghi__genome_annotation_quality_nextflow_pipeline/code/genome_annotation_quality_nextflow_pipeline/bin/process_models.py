import os
import sys
import gemmi
import numpy as np
import pandas as pd
import pickle
from multiprocessing import Pool, cpu_count
from pathlib import Path

OVERLAP_SIZE = 100
OVERLAP_PLDTT_DIFF_WARN = 20.0


def find_structure_file(predictions_dir):
    for f in os.listdir(predictions_dir):
        if f.endswith(".cif") or f.endswith(".pdb"):
            return os.path.join(predictions_dir, f)
    return None


def read_residue_plddt_ordered(file_path):
    residues = []

    if file_path.endswith(".cif"):
        cif_doc = gemmi.cif.read_file(file_path)
        block = cif_doc.sole_block()
        atom_site = block.get_mmcif_category("_atom_site")

        current_res = None
        bvals = []

        for asym, seq, b in zip(
            atom_site["label_asym_id"],
            atom_site["label_seq_id"],
            atom_site["B_iso_or_equiv"]
        ):
            res_id = (asym, int(seq))

            if current_res is None:
                current_res = res_id

            if res_id != current_res:
                residues.append(float(np.mean(bvals)))
                bvals = []
                current_res = res_id

            bvals.append(float(b))

        if bvals:
            residues.append(float(np.mean(bvals)))

    else:
        structure = gemmi.read_structure(file_path)
        for model in structure:
            for chain in model:
                for res in chain:
                    b_list = [atom.b_iso for atom in res]
                    if b_list:
                        residues.append(float(np.mean(b_list)))

    return residues


def process_single_protein(protein_dir):
    protein_id = os.path.basename(protein_dir)
    predictions_dir = os.path.join(protein_dir, "seed_101", "predictions")

    if not os.path.isdir(predictions_dir):
        return protein_id, None

    file_path = find_structure_file(predictions_dir)
    if not file_path:
        return protein_id, None

    residues = read_residue_plddt_ordered(file_path)
    if not residues:
        return protein_id, None

    return protein_id, float(np.mean(residues))


def collect_mean_plddt(dataset_path):
    dataset_name = os.path.basename(dataset_path)
    print(f"[START] Processing dataset: {dataset_name}", flush=True)

    protein_dirs = [
        os.path.join(dataset_path, d)
        for d in os.listdir(dataset_path)
        if os.path.isdir(os.path.join(dataset_path, d))
    ]

    print(f"[INFO] {dataset_name}: found {len(protein_dirs)} proteins", flush=True)

    n_workers = min(cpu_count(), 8)
    with Pool(n_workers) as pool:
        results = pool.map(process_single_protein, protein_dirs)

    plddt_dict = {pid: val for pid, val in results if val is not None}

    print(
        f"[DONE] Finished dataset: {dataset_name} "
        f"({len(plddt_dict)} proteins)\n",
        flush=True
    )

    return dataset_name, plddt_dict


def parse_args():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_models.py <dataset_all_one_dir> [output.pkl]")

    dataset_dir = Path(sys.argv[1])
    if not dataset_dir.exists():
        sys.exit(f"ERROR: dataset directory not found: {dataset_dir}")

    if len(sys.argv) >= 3:
        output_pkl = Path(sys.argv[2])
    else:
        dataset_name = dataset_dir.name
        output_pkl = Path(f"plddt_all_values_{dataset_name}.pkl")

    return dataset_dir, output_pkl


if __name__ == "__main__":
    dataset_dir, output_pkl = parse_args()
    dataset_name, plddt_dict = collect_mean_plddt(dataset_dir)

    csv_file = output_pkl.with_suffix(".csv")

    with open(output_pkl, "wb") as f:
        pickle.dump({dataset_name: plddt_dict}, f)

    rows = [
        {"Dataset": dataset_name, "Protein_ID": pid, "Mean_pLDDT": val}
        for pid, val in plddt_dict.items()
    ]

    df = pd.DataFrame(rows)
    df.to_csv(csv_file, index=False)

    print(f"Saved: {output_pkl}")
    print(f"Saved: {csv_file}")

