from __future__ import annotations
from argparse import Namespace, ArgumentParser
import sys
from ..services.batch_service import run_batch_service


def add_batch_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-i", "--cleaned_input_path",required=True,help="Path to input cleaned CIF/PDB file. The file needs to already be cleaned.")
    parser.add_argument("-m", "--input_msa",required=False,default=None,help="Optional path to input MSA file (STO/aligned FASTA/A3M format). If omitted, conservation analysis and HMM output are skipped.")
    parser.add_argument("-s","--substrate_names",required=False,default=None,help="Optional substrate names or SMILES strings. Multiple substrates should be separated by ';'. If not provided, substrate, dock, and protein-substrate interaction steps will be skipped.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to output directory. If --save_extra_outputs is False, only final integrate JSON files and log.txt will be kept here.")
    parser.add_argument("--save_extra_outputs", dest="save_extra_outputs", action="store_true",help="Enable keeping extra output files such as HMM, substrate SDFs, docked SDFs, and complex CIF files (default: Disabled).")
    parser.set_defaults(save_extra_outputs=False)
    parser.add_argument("--hydrocluster_cutoff",type=float,default=10.0,help="Minimum contact area cutoff for hydrophobic cluster residue-residue connection (default: 10.0).")
    parser.add_argument("--no_minimize_energy",action="store_false",dest="minimize_energy",help="Disable performing an energy minimization before energy evaluation (default: enabled).")
    parser.set_defaults(minimize_energy=True)
    parser.add_argument("--energy_minimization_iteration",type=int,default=100,help="Maximum number of iterations for energy minimization (default: 100).")
    parser.add_argument("--flexibility_method",type=str,choices=["ANM", "GNM"],default="ANM",help="Method for RMSF calculation: ANM or GNM (default: ANM).")
    parser.add_argument("--flexibility_cutoff",type=float,default=15.0,help="Distance cutoff used to determine the residue connection in ProDy (default: 15.0).")
    parser.add_argument("--flexibility_n_modes",type=int,default=20,help="Number of low-frequency normal modes used for RMSF calculation (default: 20).")
    parser.add_argument("--disorder_window_size",type=int,default=11,help="Sliding window size for FoldIndex-like disorder score calculation (default: 11).")
    parser.add_argument("--disorder_min_region_length",type=int,default=5,help="Minimum number of consecutive residues required to define a disordered region (default: 5).")
    parser.add_argument("--embedding_model_name",type=str,choices=["esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D"],default="esm2_t6_8M_UR50D",help="Model for embedding generation: esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, esm2_t30_150M_UR50D.")
    parser.add_argument("--pocket_min_rad",type=float,default=1.8,help="Minimum probe radius used by PyVOL for cavity detection (default: 1.8).")
    parser.add_argument("--pocket_max_rad",type=float,default=6.2,help="Maximum probe radius used by PyVOL for cavity detection (default: 6.2).")
    parser.add_argument("--pocket_min_volume",type=int,default=50,help="Minimum pocket volume threshold (default: 50).")
    parser.add_argument("--substrate_max_synonyms",type=int,default=20,help="Maximum number of substrate synonyms retried when fetching SMILES from a substrate name (default: 20).")
    parser.add_argument("--substrate_fp_radius",type=int,default=2,help="Radius used for Morgan fingerprint generation (default: 2).")
    parser.add_argument("--substrate_n_bits",type=int,default=512,help="Bit size of the Morgan fingerprint vector (default: 512).")
    parser.add_argument("--substrate_num_confs",type=int,default=5,help="Maximum number of 3D structures to generate for each substrate (default: 5).")
    parser.add_argument("--substrate_prune_rms",type=float,default=0.5,help="RMS threshold used to prune highly similar conformers during 3D conformer generation (default: 0.5).")
    parser.add_argument("--dock_max_attempt_num",type=int,default=20,help="Maximum number of docking attempts (default: 20).")
    parser.add_argument("--dock_no_early_stop", action="store_false", dest="dock_early_stop",help="Disable stopping immediately after the first successful docking result (default: enabled).")
    parser.set_defaults(dock_early_stop=True)
    parser.add_argument("--dock_exhaustiveness",type=int,default=8,help="Exhaustiveness of AutoDock Vina search (default: 8).")
    parser.add_argument("--dock_cpu",type=int,default=0,help="Number of CPUs used by AutoDock Vina (default: 0).")
    parser.add_argument("--dock_catalytic_residue",type=int,default=None,help="Cleaned protein residue index used as the docking box center. Requires '--dock_box_size'. When this parameter is provided, PyVOL pocket detection is not used.")
    parser.add_argument("--dock_catalytic_site_coord",default=None,help="Catalytic site center coordinate separated by ',' as 'x,y,z'. Requires '--dock_box_size'. When this parameter is provided, PyVOL pocket detection is not used.")
    parser.add_argument("--dock_box_size",default=None,help="Docking box size separated by ',' as 'x,y,z'. Required when '--dock_catalytic_residue' or '--dock_catalytic_site_coord' is provided.")
    parser.add_argument("--hbond_bonded_h_min_distance",type=float,default=0.8,help="Minimum bonded heavy atom-hydrogen distance used for hydrogen bond donor detection (default: 0.8).")
    parser.add_argument("--hbond_bonded_h_max_distance",type=float,default=1.3,help="Maximum bonded heavy atom-hydrogen distance used for hydrogen bond donor detection (default: 1.3).")
    parser.add_argument("--hbond_da_max_distance",type=float,default=3.9,help="Maximum donor-acceptor distance cutoff for hydrogen bond detection (default: 3.9).")
    parser.add_argument("--hbond_ha_max_distance",type=float,default=2.5,help="Maximum hydrogen-acceptor distance cutoff for hydrogen bond detection (default: 2.5).")
    parser.add_argument("--hbond_angle",type=float,default=90.0,help="Minimum donor-hydrogen-acceptor angle cutoff for hydrogen bond detection (default: 90.0).")
    parser.add_argument("--ionic_distance_cutoff",type=float,default=4.0,help="Maximum distance cutoff for ionic bond detection (default: 4.0).")
    parser.add_argument("--vdw_mu",type=float,default=0.01,help="Mu parameter used in van der Waals interaction detection (default: 0.01).")
    parser.add_argument("--ppstack_center_distance_cutoff",type=float,default=6.5,help="Maximum ring-center distance cutoff for pi-pi stacking detection (default: 6.5).")
    parser.add_argument("--pication_distance_cutoff",type=float,default=5.0,help="Maximum ring-cation distance cutoff for pi-cation interaction detection (default: 5.0).")
    parser.add_argument("--pication_angle_cutoff",type=float,default=45.0,help="Maximum angle cutoff for pi-cation interaction detection (default: 45.0).")
    parser.add_argument("--ssbond_max_distance",type=float,default=2.5,help="Maximum sulfur-sulfur distance cutoff for disulfide bond detection (default: 2.5).")


    parser.set_defaults(func=run_batch)


def run_batch(args: Namespace) -> None:
    success = run_batch_service(
        cleaned_input_path=args.cleaned_input_path,
        input_msa=args.input_msa,
        substrate_names=args.substrate_names,
        output_dir=args.output_dir,
        save_extra_outputs=args.save_extra_outputs,
        cutoff_area=args.hydrocluster_cutoff,
        minimize_energy=args.minimize_energy,
        minimization_iteration=args.energy_minimization_iteration,
        flexibility_cutoff=args.flexibility_cutoff,
        n_modes=args.flexibility_n_modes,
        flexibility_method=args.flexibility_method,
        window_size=args.disorder_window_size,
        min_region_length=args.disorder_min_region_length,
        embedding_model_name=args.embedding_model_name,
        pocket_min_rad=args.pocket_min_rad,
        pocket_max_rad=args.pocket_max_rad,
        pocket_min_volume=args.pocket_min_volume,
        max_synonyms=args.substrate_max_synonyms,
        fp_radius=args.substrate_fp_radius,
        n_bits=args.substrate_n_bits,
        num_confs=args.substrate_num_confs,
        prune_rms=args.substrate_prune_rms,
        max_docking_attempt_num=args.dock_max_attempt_num,
        early_stop=args.dock_early_stop,
        exhaustiveness=args.dock_exhaustiveness,
        cpu=args.dock_cpu,
        dock_min_rad=args.pocket_min_rad,
        dock_max_rad=args.pocket_max_rad,
        dock_min_volume=args.pocket_min_volume,
        dock_catalytic_residue=args.dock_catalytic_residue,
        dock_catalytic_site_coord=args.dock_catalytic_site_coord,
        dock_box_size=args.dock_box_size,
        bonded_h_min_distance_A=args.hbond_bonded_h_min_distance,
        bonded_h_max_distance_A=args.hbond_bonded_h_max_distance,
        da_max_distance_A=args.hbond_da_max_distance,
        ha_max_distance_A=args.hbond_ha_max_distance,
        dha_min_angle_deg=args.hbond_angle,
        ionic_distance_cutoff_A=args.ionic_distance_cutoff,
        mu=args.vdw_mu,
        ring_center_distance_cutoff_A=args.ppstack_center_distance_cutoff,
        ring_cation_distance_cutoff_A=args.pication_distance_cutoff,
        ring_cation_angle_cutoff_deg=args.pication_angle_cutoff,
        ss_max_distance_A=args.ssbond_max_distance,
    )
    if not success:
        sys.exit(1)
