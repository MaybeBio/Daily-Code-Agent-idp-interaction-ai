from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
from typing import List

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists, get_stem, check_filename_length
from ..utils.batch_utils import validate_batch_parameter_ranges, save_batch_integrate_outputs
from ..algorithms.batch_algorithms import run_batch_workflow
from ..utils.common_utils import get_optimized_filename


def _parse_float_triplet(value: str, parameter_name: str, logger: Logger) -> List[float] | None:
    if not isinstance(value, str) or not value.strip():
        logger.print(f"[ERROR] {parameter_name} is empty.")
        return None

    part_list = [x.strip() for x in value.split(",")]
    if len(part_list) != 3 or any(not x for x in part_list):
        logger.print(f"[ERROR] {parameter_name} must contain exactly 3 values separated by ','.")
        return None

    try:
        return [float(x) for x in part_list]
    except Exception as e:
        logger.print(f"[ERROR] {parameter_name} must contain numeric values. Reason: {e}")
        return None


def run_batch_service(
    cleaned_input_path: str | Path,
    input_msa: str | Path | None,
    substrate_names: str | None,
    output_dir: str | Path,
    save_extra_outputs: bool = False,
    cutoff_area: float = 10.0,
    minimize_energy: bool = True,
    minimization_iteration: int = 100,
    energy_force_field_file: str = "charmm36.xml",
    flexibility_cutoff: float = 15.0,
    n_modes: int = 20,
    flexibility_method: str = "ANM",
    window_size: int = 11,
    min_region_length: int = 5,
    embedding_model_name: str = "esm2_t6_8M_UR50D",
    pocket_min_rad: float = 1.8,
    pocket_max_rad: float = 6.2,
    pocket_min_volume: int = 50,
    max_synonyms: int = 20,
    fp_radius: int = 2,
    n_bits: int = 512,
    num_confs: int = 5,
    prune_rms: float = 0.5,
    max_docking_attempt_num: int = 20,
    early_stop: bool = True,
    exhaustiveness: int = 8,
    cpu: int = 0,
    dock_min_rad: float = 1.8,
    dock_max_rad: float = 6.2,
    dock_min_volume: int = 50,
    dock_catalytic_residue: int | None = None,
    dock_catalytic_site_coord: str | None = None,
    dock_box_size: str | None = None,
    bonded_h_min_distance_A: float = 0.8,
    bonded_h_max_distance_A: float = 1.3,
    da_max_distance_A: float = 3.9,
    ha_max_distance_A: float = 2.5,
    dha_min_angle_deg: float = 90.0,
    ionic_distance_cutoff_A: float = 4.0,
    mu: float = 0.01,
    ring_center_distance_cutoff_A: float = 6.5,
    ring_cation_distance_cutoff_A: float = 5.0,
    ring_cation_angle_cutoff_deg: float = 45.0,
    ss_max_distance_A: float = 2.5,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> bool:
    cleaned_input_path = Path(cleaned_input_path)
    has_msa = input_msa is not None and str(input_msa).strip() != ""
    input_msa_path = Path(input_msa) if has_msa else None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir_ctx: TemporaryDirectory | None = None

    try:
        if not save_extra_outputs:
            tmp_dir_ctx = TemporaryDirectory()
            working_output_dir = Path(tmp_dir_ctx.name)
        else:
            working_output_dir = output_dir

        logger = Logger(working_output_dir)

        has_substrate = isinstance(substrate_names, str) and substrate_names.strip() != ""

        logger.print(f"[INFO] Batch processing started: cleaned_input_path={cleaned_input_path}, input_msa={input_msa_path}, substrate_names={substrate_names}")
        if not has_msa:
            logger.print("[INFO] No MSA input detected. Conservation analysis and HMM output will be skipped.")
        if has_substrate:
            logger.print("[INFO] Substrate input detected. Full batch workflow will be executed.")

        dock_catalytic_site_coord_text = str(dock_catalytic_site_coord).strip() if dock_catalytic_site_coord is not None else ""
        dock_box_size_text = str(dock_box_size).strip() if dock_box_size is not None else ""
        use_manual_docking_box = dock_catalytic_residue is not None or bool(dock_catalytic_site_coord_text)

        if dock_catalytic_residue is not None and dock_catalytic_site_coord_text:
            logger.print("[ERROR] --dock_catalytic_residue and --dock_catalytic_site_coord cannot be used together.")
            return False

        if use_manual_docking_box and not dock_box_size_text:
            logger.print("[ERROR] --dock_box_size is required when --dock_catalytic_residue or --dock_catalytic_site_coord is provided.")
            return False

        dock_catalytic_site_coord_list: List[float] | None = None
        dock_box_size_list: List[float] | None = None

        if dock_catalytic_residue is not None and dock_catalytic_residue <= 0:
            logger.print("[ERROR] --dock_catalytic_residue must be a positive integer aa_id.")
            return False

        if dock_catalytic_site_coord_text:
            dock_catalytic_site_coord_list = _parse_float_triplet(dock_catalytic_site_coord_text, "--dock_catalytic_site_coord", logger)
            if dock_catalytic_site_coord_list is None:
                return False

        if use_manual_docking_box:
            dock_box_size_list = _parse_float_triplet(dock_box_size_text, "--dock_box_size", logger)
            if dock_box_size_list is None:
                return False
            if any(x <= 0.0 for x in dock_box_size_list):
                logger.print("[ERROR] --dock_box_size values must be positive.")
                return False

        if not file_exists(cleaned_input_path):
            logger.print(f"[ERROR] Input cleaned protein file not found: {cleaned_input_path}")
            return False

        if has_msa:
            if input_msa_path is None or not file_exists(input_msa_path):
                logger.print(f"[ERROR] Input MSA file not found: {input_msa_path}")
                return False

        protein_name = get_stem(cleaned_input_path)
        if not check_filename_length(protein_name, logger):
            return False
        logger.print(f"[INFO] Protein name resolved: {protein_name}")

        msa_name = None
        if has_msa:
            msa_name = get_stem(input_msa_path)
            if not check_filename_length(msa_name, logger):
                return False
            logger.print(f"[INFO] MSA file name resolved: {msa_name}")

        if not validate_batch_parameter_ranges(
            logger=logger,
            cutoff_area=cutoff_area,
            minimize_energy=minimize_energy,
            minimization_iteration=minimization_iteration,
            flexibility_cutoff=flexibility_cutoff,
            n_modes=n_modes,
            window_size=window_size,
            min_region_length=min_region_length,
            pocket_min_rad=pocket_min_rad,
            pocket_max_rad=pocket_max_rad,
            pocket_min_volume=pocket_min_volume,
            max_synonyms=max_synonyms,
            fp_radius=fp_radius,
            n_bits=n_bits,
            num_confs=num_confs,
            prune_rms=prune_rms,
            max_docking_attempt_num=max_docking_attempt_num,
            exhaustiveness=exhaustiveness,
            dock_min_rad=dock_min_rad,
            dock_max_rad=dock_max_rad,
            dock_min_volume=dock_min_volume,
            use_manual_docking_box=use_manual_docking_box,
            bonded_h_min_distance_A=bonded_h_min_distance_A,
            bonded_h_max_distance_A=bonded_h_max_distance_A,
            da_max_distance_A=da_max_distance_A,
            ha_max_distance_A=ha_max_distance_A,
            dha_min_angle_deg=dha_min_angle_deg,
            ionic_distance_cutoff_A=ionic_distance_cutoff_A,
            mu=mu,
            ring_center_distance_cutoff_A=ring_center_distance_cutoff_A,
            ring_cation_distance_cutoff_A=ring_cation_distance_cutoff_A,
            ring_cation_angle_cutoff_deg=ring_cation_angle_cutoff_deg,
            ss_max_distance_A=ss_max_distance_A,
            docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
            min_residue_index_gap=min_residue_index_gap,
        ):
            return False

        batch_result = run_batch_workflow(
            cleaned_input_path=cleaned_input_path,
            input_msa=input_msa_path,
            substrate_names=substrate_names,
            protein_name=protein_name,
            msa_name=msa_name,
            output_dir=working_output_dir,
            logger=logger,
            cutoff_area=cutoff_area,
            minimize_energy=minimize_energy,
            minimization_iteration=minimization_iteration,
            energy_force_field_file=energy_force_field_file,
            flexibility_cutoff=flexibility_cutoff,
            n_modes=n_modes,
            flexibility_method=flexibility_method,
            window_size=window_size,
            min_region_length=min_region_length,
            embedding_model_name=embedding_model_name,
            pocket_min_rad=pocket_min_rad,
            pocket_max_rad=pocket_max_rad,
            pocket_min_volume=pocket_min_volume,
            max_synonyms=max_synonyms,
            fp_radius=fp_radius,
            n_bits=n_bits,
            num_confs=num_confs,
            prune_rms=prune_rms,
            max_docking_attempt_num=max_docking_attempt_num,
            early_stop=early_stop,
            exhaustiveness=exhaustiveness,
            cpu=cpu,
            dock_min_rad=dock_min_rad,
            dock_max_rad=dock_max_rad,
            dock_min_volume=dock_min_volume,
            dock_catalytic_residue=dock_catalytic_residue,
            dock_catalytic_site_coord_list=dock_catalytic_site_coord_list,
            dock_box_size_list=dock_box_size_list,
            bonded_h_min_distance_A=bonded_h_min_distance_A,
            bonded_h_max_distance_A=bonded_h_max_distance_A,
            da_max_distance_A=da_max_distance_A,
            ha_max_distance_A=ha_max_distance_A,
            dha_min_angle_deg=dha_min_angle_deg,
            ionic_distance_cutoff_A=ionic_distance_cutoff_A,
            mu=mu,
            ring_center_distance_cutoff_A=ring_center_distance_cutoff_A,
            ring_cation_distance_cutoff_A=ring_cation_distance_cutoff_A,
            ring_cation_angle_cutoff_deg=ring_cation_angle_cutoff_deg,
            ss_max_distance_A=ss_max_distance_A,
            docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
            min_residue_index_gap=min_residue_index_gap,
        )
        if batch_result is None:
            return False

        integrate_report = batch_result.get("integrate_report")
        if not isinstance(integrate_report, dict):
            logger.print("[ERROR] Missing integrate_report in batch result.")
            return False

        if not save_batch_integrate_outputs(
            integrate_report=integrate_report,
            output_dir=working_output_dir,
            protein_name=protein_name,
            logger=logger,
        ):
            return False

        report_filename = get_optimized_filename(f"integrate_report_{protein_name}.json")
        nodes_filename = get_optimized_filename(f"integrate_nodes_{protein_name}.json")
        edges_filename = get_optimized_filename(f"integrate_edges_{protein_name}.json")
        log_filename = "log.txt"

        logger.print("[INFO] Batch processing finished")

        if not save_extra_outputs:
            shutil.copy2(working_output_dir / report_filename, output_dir / report_filename)
            shutil.copy2(working_output_dir / nodes_filename, output_dir / nodes_filename)
            shutil.copy2(working_output_dir / edges_filename, output_dir / edges_filename)

            log_path = working_output_dir / log_filename
            if log_path.exists():
                shutil.copy2(log_path, output_dir / log_filename)

        return True

    finally:
        if tmp_dir_ctx is not None:
            tmp_dir_ctx.cleanup()
