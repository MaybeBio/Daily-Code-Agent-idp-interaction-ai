####################################
############## ColabDesign functions
####################################
### Import dependencies
import os, re, shutil, math, pickle
import matplotlib.pyplot as plt
import numpy as np
import jax
import jax.numpy as jnp
from scipy.special import softmax
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.loss import get_ptm, get_pae, mask_loss, get_dgram_bins, _get_con_loss
from colabdesign.shared.utils import copy_dict
from .biopython_utils import hotspot_residues, calculate_clash_score, calc_ss_percentage, calculate_percentages, interface_residue_indices, calculate_clashes_after_superpostion_on_full_target, superimpose_on_full_target_and_save_pdb, binder_crop_terminus_engagement, binder_terminus_clearance
from .pyrosetta_utils import pr_relax, align_pdbs
from .generic_utils import update_failures
from .ipsae import compute_ipsae

# hallucinate a binder
def binder_hallucination(design_name, starting_pdb, chain, target_hotspot_residues, length, seed, helicity_value, design_models, advanced_settings, design_paths, failure_csv, full_target_structure=None, full_target_pae=None, crop_structure_indices=None):
    model_pdb_path = os.path.join(design_paths["Trajectory"], design_name+".pdb")

    # clear GPU memory for new trajectory
    clear_mem()

    # initialise binder hallucination model
    af_model = mk_afdesign_model(protocol="binder", debug=False, data_dir=advanced_settings["af_params_dir"], 
                                use_multimer=advanced_settings["use_multimer_design"], num_recycles=advanced_settings["num_recycles_design"],
                                best_metric='loss')

    # sanity check for hotspots
    if target_hotspot_residues == "":
        target_hotspot_residues = None

    af_model.prep_inputs(pdb_filename=starting_pdb, chain=chain, binder_len=length, hotspot=target_hotspot_residues, seed=seed, rm_aa=advanced_settings["omit_AAs"],
                        rm_target=advanced_settings.get("rm_template_target_design", False),
                        rm_target_seq=advanced_settings["rm_template_seq_design"], rm_target_sc=advanced_settings["rm_template_sc_design"])

    ### Update weights based on specified settings
    af_model.opt["weights"].update({"pae":advanced_settings["weights_pae_intra"],
                                    "plddt":advanced_settings["weights_plddt"],
                                    "i_pae":advanced_settings["weights_pae_inter"],
                                    "con":advanced_settings["weights_con_intra"],
                                    "i_con":advanced_settings["weights_con_inter"],
                                    })

    # redefine intramolecular contacts (con) and intermolecular contacts (i_con) definitions
    af_model.opt["con"].update({"num":advanced_settings["intra_contact_number"],"cutoff":advanced_settings["intra_contact_distance"],"binary":False,"seqsep":9})
    af_model.opt["i_con"].update({"num":advanced_settings["inter_contact_number"],"cutoff":advanced_settings["inter_contact_distance"],"binary":False})
        

    ### additional loss functions
    if advanced_settings["use_rg_loss"]:
        # radius of gyration loss
        add_rg_loss(af_model, advanced_settings["weights_rg"])

    if advanced_settings["use_i_ptm_loss"]:
        # interface pTM loss
        add_i_ptm_loss(af_model, advanced_settings["weights_iptm"])

    if advanced_settings["use_termini_distance_loss"]:
        # termini distance loss
        add_termini_distance_loss(af_model, advanced_settings["weights_termini_loss"])

    if advanced_settings.get("use_interface_plddt_target_loss", False):
        # IDP-specific: pLDDT only on contacted target residues (additive)
        add_interface_plddt_target_loss(af_model, advanced_settings.get("weights_plddt_target_iface", 0.1))

    if advanced_settings.get("use_interface_pae_loss", False):
        add_interface_pae_loss(af_model, advanced_settings.get("weights_pae_iface", 0.1))

    if advanced_settings.get("use_actifptm_loss", False):
        # actifpTM: contact-probability-weighted ipTM (weight 0 => metric only)
        add_actifptm_loss(af_model,
                          advanced_settings.get("weights_actifptm", 0.0),
                          advanced_settings.get("actifptm_contact_dist", 8.0))

    if advanced_settings.get("use_actif_i_ptm_energy_loss", False):
        # hybrid: dense pTMEnergy restricted to actually-contacting pairs
        add_actif_i_ptm_energy_loss(af_model,
                                    advanced_settings.get("weights_actif_iptmenergy_loss", 1.0),
                                    advanced_settings.get("iptmenergy_lambda", 1.0),
                                    advanced_settings.get("actifptm_contact_dist", 8.0))

    if advanced_settings.get("use_i_ptm_energy_loss", False):
        # BECraft: interface pTMEnergy loss, a REPLACEMENT for i_ptm (not additive).
        # Both can be registered at once, which optimises two interface objectives on
        # different scales (i_ptm in [0,1], the energy unbounded and negative), so warn
        # loudly rather than silently doing something the config author likely did not mean.
        if advanced_settings["use_i_ptm_loss"] and advanced_settings["weights_iptm"] != 0:
            print("WARNING: both use_i_ptm_loss and use_i_ptm_energy_loss are enabled "
                  f"(weights_iptm={advanced_settings['weights_iptm']}, "
                  f"weights_iptmenergy_loss={advanced_settings.get('weights_iptmenergy_loss', 1.0)}). "
                  "pTMEnergy is intended to REPLACE ipTM; set use_i_ptm_loss false and "
                  "weights_iptm 0 unless you deliberately want both.")
        add_i_ptm_energy_loss(af_model,
                              advanced_settings.get("weights_iptmenergy_loss", 1.0),
                              advanced_settings.get("iptmenergy_lambda", 1.0))

    # add the helicity loss
    add_helix_loss(af_model, helicity_value)

    # calculate the number of mutations to do based on the length of the protein
    greedy_tries = math.ceil(length * (advanced_settings["greedy_percentage"] / 100))

    ### start design algorithm based on selection
    if advanced_settings["design_algorithm"] == '2stage':
        # uses gradient descend to get a PSSM profile and then uses PSSM to bias the sampling of random mutations to decrease loss
        af_model.design_pssm_semigreedy(soft_iters=advanced_settings["soft_iterations"], hard_iters=advanced_settings["greedy_iterations"], tries=greedy_tries, models=design_models, 
                                        num_models=1, sample_models=advanced_settings["sample_models"], ramp_models=False, save_best=True)

    elif advanced_settings["design_algorithm"] == '3stage':
        # 3 stage design using logits, softmax, and one hot encoding
        af_model.design_3stage(soft_iters=advanced_settings["soft_iterations"], temp_iters=advanced_settings["temporary_iterations"], hard_iters=advanced_settings["hard_iterations"], 
                                num_models=1, models=design_models, sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == 'greedy':
        # design by using random mutations that decrease loss
        af_model.design_semigreedy(advanced_settings["greedy_iterations"], tries=greedy_tries, num_models=1, models=design_models,
                                sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == 'mcmc':
        # design by using random mutations that decrease loss
        half_life = round(advanced_settings["greedy_iterations"] / 5, 0)
        t_mcmc = 0.01
        af_model._design_mcmc(advanced_settings["greedy_iterations"], half_life=half_life, T_init=t_mcmc, mutation_rate=greedy_tries, num_models=1, models=design_models,
                                sample_models=advanced_settings["sample_models"], save_best=True)

    elif advanced_settings["design_algorithm"] == '4stage':
        # initial logits to prescreen trajectory
        print("Stage 1: Test Logits")
        af_model.design_logits(iters=50, e_soft=0.9, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True)

        # determine pLDDT of best iteration according to lowest 'loss' value
        initial_plddt = get_best_plddt(af_model, length)
        
        # if best iteration has high enough confidence then continue
        if initial_plddt > 0.65:
            print("Initial trajectory pLDDT good, continuing: "+str(initial_plddt))
            if advanced_settings["optimise_beta"]:
                # temporarily dump model to assess secondary structure
                af_model.save_pdb(model_pdb_path)
                _, beta, *_ = calc_ss_percentage(model_pdb_path, advanced_settings, 'B')
                os.remove(model_pdb_path)

                # if beta sheeted trajectory is detected then choose to optimise
                if float(beta) > 15:
                    advanced_settings["soft_iterations"] = advanced_settings["soft_iterations"] + advanced_settings["optimise_beta_extra_soft"]
                    advanced_settings["temporary_iterations"] = advanced_settings["temporary_iterations"] + advanced_settings["optimise_beta_extra_temp"]
                    af_model.set_opt(num_recycles=advanced_settings["optimise_beta_recycles_design"])
                    print("Beta sheeted trajectory detected, optimising settings")

            # how many logit iterations left
            logits_iter = advanced_settings["soft_iterations"] - 50
            if logits_iter > 0:
                print("Stage 1: Additional Logits Optimisation")
                af_model.clear_best()
                af_model.design_logits(iters=logits_iter, e_soft=1, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"],
                                    ramp_recycles=False, save_best=True)
                af_model._tmp["seq_logits"] = af_model.aux["seq"]["logits"]
                logit_plddt = get_best_plddt(af_model, length)
                print("Optimised logit trajectory pLDDT: "+str(logit_plddt))
            else:
                logit_plddt = initial_plddt

            # perform softmax trajectory design
            if advanced_settings["temporary_iterations"] > 0:
                print("Stage 2: Softmax Optimisation")
                af_model.clear_best()
                af_model.design_soft(advanced_settings["temporary_iterations"], e_temp=1e-2, models=design_models, num_models=1,
                                    sample_models=advanced_settings["sample_models"], ramp_recycles=False, save_best=True)
                softmax_plddt = get_best_plddt(af_model, length)
            else:
                softmax_plddt = logit_plddt

            # perform one hot encoding
            if softmax_plddt > 0.65:
                print("Softmax trajectory pLDDT good, continuing: "+str(softmax_plddt))

                onehot_plddt = softmax_plddt

                if advanced_settings["hard_iterations"] > 0:
                    af_model.clear_best()
                    print("Stage 3: One-hot Optimisation")
                    af_model.design_hard(advanced_settings["hard_iterations"], temp=1e-2, models=design_models, num_models=1,
                                    sample_models=advanced_settings["sample_models"], dropout=False, ramp_recycles=False, save_best=True)
                    onehot_plddt = get_best_plddt(af_model, length)

                if onehot_plddt > 0.65:
                    # perform greedy mutation optimisation
                    print("One-hot trajectory pLDDT good, continuing: "+str(onehot_plddt))
                    if advanced_settings["greedy_iterations"] > 0:
                        print("Stage 4: PSSM Semigreedy Optimisation")
                        af_model.design_pssm_semigreedy(soft_iters=0, hard_iters=advanced_settings["greedy_iterations"], tries=greedy_tries, models=design_models, 
                                                        num_models=1, sample_models=advanced_settings["sample_models"], ramp_models=False, save_best=True)

                else:
                    update_failures(failure_csv, 'Trajectory_one-hot_pLDDT')
                    print("One-hot trajectory pLDDT too low to continue: "+str(onehot_plddt))

            else:
                update_failures(failure_csv, 'Trajectory_softmax_pLDDT')
                print("Softmax trajectory pLDDT too low to continue: "+str(softmax_plddt))

        else:
            update_failures(failure_csv, 'Trajectory_logits_pLDDT')
            print("Initial trajectory pLDDT too low to continue: "+str(initial_plddt))

    else:
        print("ERROR: No valid design model selected")
        exit()
        return

    ### save trajectory PDB
    final_plddt = get_best_plddt(af_model, length)
    af_model.save_pdb(model_pdb_path)
    af_model.aux["log"]["terminate"] = ""

    # let's check whether the trajectory is worth optimising by checking confidence, clashes, and contacts
    # check clashes
    #clash_interface = calculate_clash_score(model_pdb_path, 2.4)
    ca_clashes = calculate_clash_score(model_pdb_path, 2.5, only_ca=True)

    #if clash_interface > 25 or ca_clashes > 0:
    if ca_clashes > 0:
        af_model.aux["log"]["terminate"] = "Clashing"
        update_failures(failure_csv, 'Trajectory_Clashes')
        print("Severe clashes detected, skipping analysis and MPNN optimisation")
        print("")
    else:
        # check if low quality prediction
        if final_plddt < 0.7:
            af_model.aux["log"]["terminate"] = "LowConfidence"
            update_failures(failure_csv, 'Trajectory_final_pLDDT')
            print("Trajectory starting confidence low, skipping analysis and MPNN optimisation")
            print("")
        else:
            # does it have enough contacts to consider?
            binder_contacts = hotspot_residues(model_pdb_path)
            binder_contacts_n = len(binder_contacts.items())

            # if less than 3 contacts then protein is floating above and is not binder
            if binder_contacts_n < 3:
                af_model.aux["log"]["terminate"] = "LowConfidence"
                update_failures(failure_csv, 'Trajectory_Contacts')
                print("Too few contacts at the interface, skipping analysis and MPNN optimisation")
                print("")
            else:
                terminate_reason = ""

                # full-target validation (opt-in): reject trajectories whose binder
                # confidently clashes with the full-length target outside the crop
                if full_target_structure is not None:
                    if advanced_settings.get("save_superimpose_trajectory_structure", False):
                        superimposed_filepath = os.path.join(design_paths["Trajectory/Superimposed"], design_name+".pdb")
                    else:
                        superimposed_filepath = None
                    clash_dict = calculate_clashes_after_superpostion_on_full_target(
                        model_pdb_path, full_target_structure, full_target_pae,
                        crop_structure_indices, superimposed_filepath=superimposed_filepath,
                        clash_radius=1.5, clash_plddt_threshold=70.0, clash_pae_threshold=15.0)
                    if clash_dict["clash_pass"] is False:
                        terminate_reason = "ClashingFullTarget"
                        update_failures(failure_csv, 'Trajectory_Clashes_FullTarget')
                        print("Clashes detected with full target structure, skipping analysis and MPNN optimisation")
                        print("")

                # crop-terminus engagement (opt-in): flag/reject binders that grip the
                # crop's terminal (cap) residues. Under full-length superimposition a cut
                # terminus is where the real chain continues, so binding it is a cropping
                # artifact. Metric is logged whenever measured; rejection is separate.
                if terminate_reason == "" and advanced_settings.get("measure_terminus_engagement", False):
                    term = binder_crop_terminus_engagement(
                        model_pdb_path,
                        target_crop_indices=crop_structure_indices,
                        full_target_structure=full_target_structure,
                        terminus_len=advanced_settings.get("terminus_length", 5),
                        only_cut_ends=advanced_settings.get("terminus_only_cut_ends", True))
                    af_model.aux["log"]["n_terminus_iface"] = term["n_terminus_iface"]
                    af_model.aux["log"]["terminus_iface_fraction"] = term["terminus_iface_fraction"]
                    # reject when the binder grips a crop CUT end: either an edge-dominated
                    # interface (fraction gate, cheap false-positive filter for extremes) or
                    # any absolute count on the cap (max_terminus_iface_count). The count
                    # gate matters because for IDP crops with ~25-50 interface residues even
                    # a handful of cap contacts can be a cropping artifact while sitting far
                    # below any sensible fraction. Count defaults to a very large number so
                    # the behaviour is unchanged unless the setting is set.
                    if advanced_settings.get("reject_terminus_binders", False):
                        max_frac = advanced_settings.get("max_terminus_iface_fraction", 0.5)
                        max_count = advanced_settings.get("max_terminus_iface_count", 10**9)
                        if (term["terminus_iface_fraction"] > max_frac
                                or term["n_terminus_iface"] > max_count):
                            terminate_reason = "BindsCropTerminus"
                            update_failures(failure_csv, 'Trajectory_BindsCropTerminus')
                            print(f"Binder grips crop cut end "
                                  f"({term['n_terminus_iface']} cap residues, "
                                  f"{round(term['terminus_iface_fraction']*100)}% of interface), "
                                  f"skipping analysis and MPNN optimisation")
                            print("")

                af_model.aux["log"]["terminate"] = terminate_reason
                if terminate_reason == "":
                    # phew, trajectory is okay! We can continue
                    print("Trajectory successful, final pLDDT: "+str(final_plddt))

    # move low quality prediction:
    if af_model.aux["log"]["terminate"] != "":
        # On a seed collision (same filename already in the reject bucket from an earlier leg),
        # shutil.move raises Error and crashes bindcraft -- which kills the SLURM resubmit chain
        # because rc=1 isn't in the walltime-signal set. Overwrite instead; the two files are
        # the same trajectory re-sampled.
        _dst_dir = design_paths[f"Trajectory/{af_model.aux['log']['terminate']}"]
        _dst = os.path.join(_dst_dir, os.path.basename(model_pdb_path))
        try:
            if os.path.exists(_dst):
                os.remove(_dst)
            shutil.move(model_pdb_path, _dst_dir)
        except Exception as _e:
            print(f"[trajectory move] WARNING: could not move {os.path.basename(model_pdb_path)} to "
                  f"{af_model.aux['log']['terminate']} ({type(_e).__name__}: {_e}); discarding file.")
            try:
                os.remove(model_pdb_path)
            except OSError:
                pass

    ### get the sampled sequence for plotting
    af_model.get_seqs()
    if advanced_settings["save_design_trajectory_plots"]:
        plot_trajectory(af_model, design_name, design_paths)

    ### save the hallucination trajectory animation
    if advanced_settings["save_design_animations"]:
        plots = af_model.animate(dpi=150)
        with open(os.path.join(design_paths["Trajectory/Animation"], design_name+".html"), 'w') as f:
            f.write(plots)
        plt.close('all')

    if advanced_settings["save_trajectory_pickle"]:
        with open(os.path.join(design_paths["Trajectory/Pickle"], design_name+".pickle"), 'wb') as handle:
            pickle.dump(af_model.aux['all'], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return af_model

# run prediction for binder with masked template target
def predict_binder_complex(prediction_model, binder_sequence, mpnn_design_name, target_pdb, chain, length, trajectory_pdb, prediction_models, advanced_settings, filters, design_paths, failure_csv, seed=None):
    prediction_stats = {}

    # clean sequence
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())

    # reset filtering conditionals
    pass_af2_filters = True
    filter_failures = {}

    # start prediction per AF2 model, 2 are used by default due to masked templates
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if not os.path.exists(complex_pdb):
            # predict model
            prediction_model.predict(seq=binder_sequence, models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            prediction_model.save_pdb(complex_pdb)
            prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

            # interface-only i_pAE and contact-restricted i_pTM (project convention: 4 Å all-atom)
            # Both restrict to residues in actual contact, so floppy non-contact target
            # residues (IDP tails) do not dilute the metric. Consistent with the interface loss.
            iface_iptm_val = float(prediction_metrics["i_ptm"])
            n_term_iface_val = 0
            term_frac_val = 0.0
            try:
                pae_matrix = np.asarray(prediction_model.aux["pae"])
                t_idx, b_idx, target_len, _ = interface_residue_indices(complex_pdb, binder_chain="B", cutoff=4.0)
                if t_idx and b_idx:
                    b_idx_abs = [target_len + i for i in b_idx]
                    pae_tb = pae_matrix[np.ix_(t_idx, b_idx_abs)]
                    pae_bt = pae_matrix[np.ix_(b_idx_abs, t_idx)]
                    iface_pae_val = float((pae_tb.mean() + pae_bt.mean()) / 2 / 31.0)

                    # contact-restricted i_pTM: pTM formula over interface residues only,
                    # anchor over each interface residue, summing cross-chain interface pairs
                    iface_abs = list(t_idx) + b_idx_abs
                    n_if = len(iface_abs)
                    d0 = max(1.24 * (max(n_if, 19) - 15) ** (1.0 / 3.0) - 1.8, 0.5)
                    best = 0.0
                    for i in iface_abs:
                        js = b_idx_abs if i < target_len else t_idx   # cross-chain interface residues
                        paes = pae_matrix[i, js]
                        s = float((1.0 / (1.0 + (paes / d0) ** 2)).mean())
                        if s > best:
                            best = s
                    iface_iptm_val = best

                    # crop-terminus engagement: how many target interface residues fall in
                    # the first/last K crop residues (the hydrophilic caps). Positional, so
                    # numbering-independent. Both ends counted here (the cut-vs-real-terminus
                    # refinement lives in the trajectory stage where the full structure is loaded).
                    K = min(int(advanced_settings.get("terminus_length", 5)), target_len)
                    term_pos = set(range(0, K)) | set(range(target_len - K, target_len))
                    n_term_iface_val = sum(1 for i in t_idx if i in term_pos)
                    term_frac_val = round(n_term_iface_val / len(t_idx), 3) if t_idx else 0.0
                else:
                    iface_pae_val = float(prediction_metrics["i_pae"])
            except Exception:
                iface_pae_val = float(prediction_metrics["i_pae"])
            prediction_metrics["i_pae_iface"] = iface_pae_val
            prediction_metrics["i_ptm_iface"] = iface_iptm_val
            prediction_metrics["n_terminus_iface"] = n_term_iface_val
            prediction_metrics["terminus_iface_fraction"] = term_frac_val

            # binder-terminus clearance: room at the binder's own N/C terminus for a tag
            # or fusion. Distinct from terminus_iface_fraction above, which is about the
            # TARGET crop's cut ends. Long binders wrapping a short extended crop score
            # well on every confidence metric yet bury both termini, making them untaggable.
            term_clear_val = 0.0
            if advanced_settings.get("measure_binder_terminus_clearance", False):
                try:
                    term_clear_val = binder_terminus_clearance(
                        complex_pdb,
                        iface_cutoff=advanced_settings.get("binder_terminus_iface_cutoff", 5.0),
                    )["binder_term_clearance"]
                except Exception:
                    term_clear_val = 0.0
            prediction_metrics["binder_term_clearance"] = term_clear_val

            # ipSAE (Dunbrack, bioRxiv 2025.02.10.637595). Selects interface pairs by PAE
            # alone and normalises d0 per residue, so unlike i_pTM_iface above it needs no
            # contact-distance definition and does not get diluted by the flexible parts of
            # an IDP crop. Computed in its own try block so it can never perturb the
            # existing metrics, and only gates if a threshold is set in the filters JSON.
            ipsae_val = 0.0
            try:
                ipsae_pae = np.asarray(prediction_model.aux["pae"])
                ipsae_target_len = ipsae_pae.shape[0] - length
                if ipsae_target_len > 0:
                    # AF2 binder protocol orders residues target-first, then binder.
                    ipsae_chains = np.array(["A"] * ipsae_target_len + ["B"] * length)
                    ipsae_val = float(compute_ipsae(
                        ipsae_pae, ipsae_chains, "A", "B",
                        pae_cutoff=advanced_settings.get("ipsae_pae_cutoff", 10.0),
                    )["ipsae"])
            except Exception:
                ipsae_val = 0.0
            prediction_metrics["ipsae"] = ipsae_val

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2),
                'pTM': round(prediction_metrics['ptm'], 2),
                'i_pTM': round(prediction_metrics['i_ptm'], 2),
                'pAE': round(prediction_metrics['pae'], 2),
                'i_pAE': round(prediction_metrics['i_pae'], 2),
                'i_pAE_iface': round(iface_pae_val, 2),
                'i_pTM_iface': round(iface_iptm_val, 2),
                'n_terminus_iface': n_term_iface_val,
                'terminus_iface_fraction': round(term_frac_val, 2),
                'binder_term_clearance': round(term_clear_val, 2),
                # 3 dp: ipSAE is used to rank designs, and 2 dp collapses too many ties
                'ipSAE': round(ipsae_val, 3)
            }
            prediction_stats[model_num+1] = stats

            # List of filter conditions and corresponding keys
            filter_conditions = [
                (f"{model_num+1}_pLDDT", 'plddt', '>='),
                (f"{model_num+1}_pTM", 'ptm', '>='),
                (f"{model_num+1}_i_pTM", 'i_ptm', '>='),
                (f"{model_num+1}_pAE", 'pae', '<='),
                (f"{model_num+1}_i_pAE", 'i_pae', '<='),
                (f"{model_num+1}_i_pAE_iface", 'i_pae_iface', '<='),
                (f"{model_num+1}_i_pTM_iface", 'i_ptm_iface', '>='),
                (f"{model_num+1}_ipSAE", 'ipsae', '>='),
                (f"{model_num+1}_terminus_iface_fraction", 'terminus_iface_fraction', '<='),
                (f"{model_num+1}_n_terminus_iface", 'n_terminus_iface', '<='),
                (f"{model_num+1}_binder_term_clearance", 'binder_term_clearance', '>='),
            ]

            # perform initial AF2 values filtering to determine whether to skip relaxation and interface scoring
            for filter_name, metric_key, comparison in filter_conditions:
                threshold = filters.get(filter_name, {}).get("threshold")
                if threshold is not None:
                    if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                    elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1

            if not pass_af2_filters:
                break

    # Update the CSV file with the failure counts
    if filter_failures:
        update_failures(failure_csv, filter_failures)

    # AF2 filters passed, contuing with relaxation
    for model_num in prediction_models:
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if pass_af2_filters:
            mpnn_relaxed = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{model_num+1}.pdb")
            pr_relax(complex_pdb, mpnn_relaxed)
        else:
            if os.path.exists(complex_pdb):
                os.remove(complex_pdb)

    return prediction_stats, pass_af2_filters

# run prediction for binder alone
def predict_binder_alone(prediction_model, binder_sequence, mpnn_design_name, length, trajectory_pdb, binder_chain, prediction_models, advanced_settings, design_paths, seed=None):
    binder_stats = {}

    # prepare sequence for prediction
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())
    prediction_model.set_seq(binder_sequence)

    # predict each model separately f
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        binder_alone_pdb = os.path.join(design_paths["MPNN/Binder"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if not os.path.exists(binder_alone_pdb):
            # predict model
            prediction_model.predict(models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            prediction_model.save_pdb(binder_alone_pdb)
            prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, pae

            # align binder model to trajectory binder
            try:
                align_pdbs(trajectory_pdb, binder_alone_pdb, binder_chain, "A")
            except Exception as e:
                print(f"Warning: alignment failed for {mpnn_design_name}_model{model_num+1}: {e}")
                
            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2)
            }
            binder_stats[model_num+1] = stats

    return binder_stats

# run MPNN to generate sequences for binders
def mpnn_gen_sequence(trajectory_pdb, binder_chain, trajectory_interface_residues, advanced_settings):
    # clear GPU memory
    clear_mem()

    # initialise MPNN model
    mpnn_model = mk_mpnn_model(backbone_noise=advanced_settings["backbone_noise"], model_name=advanced_settings["model_path"], weights=advanced_settings["mpnn_weights"])

    # check whether keep the interface generated by the trajectory or whether to redesign with MPNN
    design_chains = 'A,' + binder_chain

    if advanced_settings["mpnn_fix_interface"]:
        fixed_positions = 'A,' + trajectory_interface_residues
        fixed_positions = fixed_positions.rstrip(",")
        print("Fixing interface residues: "+trajectory_interface_residues)
    else:
        fixed_positions = 'A'

    # prepare inputs for MPNN
    mpnn_model.prep_inputs(pdb_filename=trajectory_pdb, chain=design_chains, fix_pos=fixed_positions, rm_aa=advanced_settings["omit_AAs"])

    # sample MPNN sequences in parallel
    mpnn_sequences = mpnn_model.sample(temperature=advanced_settings["sampling_temp"], num=1, batch=advanced_settings["num_seqs"])

    return mpnn_sequences

# Get pLDDT of best model
def get_best_plddt(af_model, length):
    return round(np.mean(af_model._tmp["best"]["aux"]["plddt"][-length:]),2)

# IDP-specific: contact-weighted pLDDT over target-side interface residues only.
# An IDP target has globally low pLDDT; pushes up
# confidence only on the target residues the binder actually contacts (Cα within 4 Å),
# so a binder that stabilises the IDP raises pLDDT locally even if the rest stays
# disordered. Does not touch binder pLDDT or disable the default plddt term.
def add_interface_plddt_target_loss(self, weight=0.1):
    target_len = self._target_len
    binder_len = self._binder_len
    bin_centers = (jnp.arange(50) + 0.5) / 50.0

    def loss_fn(inputs, outputs):
        ca = outputs["structure_module"]["final_atom_positions"][:,
                residue_constants.atom_order["CA"]]
        t_ca, b_ca = ca[:target_len], ca[-binder_len:]
        d_t = jnp.sqrt(((t_ca[:, None, :] - b_ca[None, :, :]) ** 2).sum(-1).min(1) + 1e-8)
        m_t = jax.nn.sigmoid(4.0 - d_t)   # soft contact weight per target residue

        plddt = (jax.nn.softmax(outputs["predicted_lddt"]["logits"], -1) * bin_centers).sum(-1)
        iface = (m_t * plddt[:target_len]).sum() / (m_t.sum() + 1e-8)
        return {"plddt_target_iface": 1.0 - iface}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["plddt_target_iface"] = weight


# Optimize inter-chain pAE only for interface residues (Cα within 4 Å of partner chain).
def add_interface_pae_loss(self, weight=0.1):
    target_len = self._target_len
    binder_len = self._binder_len

    def loss_fn(inputs, outputs):
        ca = outputs["structure_module"]["final_atom_positions"][:,
                residue_constants.atom_order["CA"]]
        t_ca, b_ca = ca[:target_len], ca[-binder_len:]
        d_b = jnp.sqrt(((b_ca[:, None, :] - t_ca[None, :, :]) ** 2).sum(-1).min(1) + 1e-8)
        d_t = jnp.sqrt(((t_ca[:, None, :] - b_ca[None, :, :]) ** 2).sum(-1).min(1) + 1e-8)
        m_b = jax.nn.sigmoid(4.0 - d_b)   # [binder]
        m_t = jax.nn.sigmoid(4.0 - d_t)   # [target]

        pae = get_pae(outputs)                             # [n, n] in Å
        pae_tb = pae[:target_len, target_len:]             # target frame → binder
        pae_bt = pae[target_len:, :target_len]             # binder frame → target
        mask_tb = m_t[:, None] * m_b[None, :]             # [target, binder]

        iface_pae = ((mask_tb * pae_tb).sum() + (mask_tb.T * pae_bt).sum()) \
                    / (2 * mask_tb.sum() + 1e-8)
        return {"pae_iface": iface_pae / 31.0}            # normalize to [0,1]

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["pae_iface"] = weight


# Define radius of gyration loss for colabdesign
def add_rg_loss(self, weight=0.1):
    '''add radius of gyration loss'''
    def loss_fn(inputs, outputs):
        xyz = outputs["structure_module"]
        ca = xyz["final_atom_positions"][:,residue_constants.atom_order["CA"]]
        ca = ca[-self._binder_len:]
        rg = jnp.sqrt(jnp.square(ca - ca.mean(0)).sum(-1).mean() + 1e-8)
        rg_th = 2.38 * ca.shape[0] ** 0.365

        rg = jax.nn.elu(rg - rg_th)
        return {"rg":rg}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["rg"] = weight

# Define interface pTM loss for colabdesign
def add_i_ptm_loss(self, weight=0.1):
    def loss_iptm(inputs, outputs):
        p = 1 - get_ptm(inputs, outputs, interface=True)
        i_ptm = mask_loss(p)
        return {"i_ptm": i_ptm}
    
    self._callbacks["model"]["loss"].append(loss_iptm)
    self.opt["weights"]["i_ptm"] = weight

# Define interface pTMEnergy loss for colabdesign (BECraft, arXiv:2505.21241)
# Reinterprets the pAE head as an energy-based model: instead of normalising the
# pAE logits with a softmax (as pTM does) we take -logsumexp over the TM-weighted
# logits. Dropping the partition function keeps the logit magnitudes, which gives
# a denser gradient than ipTM on hard/low-confidence interfaces.
def get_iptmenergy(inputs, outputs, tm_lambda=1.0):
    logits = outputs["predicted_aligned_error"]["logits"]   # [n, n, bins]
    breaks = outputs["predicted_aligned_error"]["breaks"]
    asym_id = inputs["asym_id"]

    # bin centres, computed the same way as colabdesign's get_pae
    step = breaks[1] - breaks[0]
    bin_centers = breaks + step / 2
    bin_centers = jnp.append(bin_centers, bin_centers[-1] + step)

    # TM-score kernel g(d_b), up-weighting small predicted errors
    num_res = jnp.maximum(logits.shape[0], 19)
    d0 = 1.24 * (num_res - 15) ** (1.0 / 3) - 1.8
    g = 1.0 / (1 + jnp.square(bin_centers) / jnp.square(d0))

    # tm_lambda scales the TM-weighting: 1 = full (default), 0 = plain -logsumexp
    weighted_logits = logits + tm_lambda * jnp.log(g)[None, None, :]
    positional_energy = -jax.nn.logsumexp(weighted_logits, axis=-1)

    # average over inter-chain pairs only
    mask = (asym_id[:, None] != asym_id[None, :]).astype(jnp.float32)
    denom = jnp.maximum(mask.sum(), 1.0)
    return (positional_energy * mask).sum() / denom


def add_i_ptm_energy_loss(self, weight=1.0, tm_lambda=1.0):
    def loss_iptmenergy(inputs, outputs):
        # already an energy (lower is better), so no "1 minus" as with i_ptm
        return {"i_ptm_energy": get_iptmenergy(inputs, outputs, tm_lambda)}

    self._callbacks["model"]["loss"].append(loss_iptmenergy)
    self.opt["weights"]["i_ptm_energy"] = weight


# ---------------------------------------------------------------------------
# actifpTM (Bioinformatics 2025, btaf107; ColabFold colabfold/alphafold/extra_ptm.py)
#
# Plain ipTM averages over EVERY inter-chain residue pair, so flexible flanking
# regions that stay unbound still drag the score down -- and two models of the same
# motif-domain interaction differing only in flank length score very differently.
# actifpTM fixes that by weighting each pair by its CONTACT PROBABILITY instead of
# uniformly, so the score reflects the residues actually participating.
# That dilution problem is exactly what IDR crops have.
# ---------------------------------------------------------------------------

def _interchain_contact_weights(inputs, outputs, contact_dist=8.0):
    """Contact-probability map from the distogram, masked to inter-chain pairs.
    Differentiable: softmax over distogram bins, same as colabdesign's get_contact_map."""
    dgram = outputs["distogram"]["logits"]
    dgram_bins = get_dgram_bins(outputs)
    cmap = (jax.nn.softmax(dgram) * (dgram_bins < contact_dist)).sum(-1)
    asym_id = inputs["asym_id"]
    pair_mask = (asym_id[:, None] != asym_id[None, :]).astype(cmap.dtype)
    return cmap * pair_mask, pair_mask


def get_actifptm(inputs, outputs, contact_dist=8.0):
    """actifpTM, matching predicted_tm_score_modified(pair_residue_weights=cmap).max()."""
    logits = outputs["predicted_aligned_error"]["logits"]
    breaks = outputs["predicted_aligned_error"]["breaks"]
    residue_weights = inputs["seq_mask"]

    step = breaks[1] - breaks[0]
    bin_centers = jnp.append(breaks + step / 2, breaks[-1] + step / 2 + step)

    # reference uses residue_weights.shape[0] (total length), not its sum
    num_res = residue_weights.shape[0]
    d0 = 1.24 * (max(num_res, 19) - 15) ** (1.0 / 3) - 1.8
    tm_per_bin = 1.0 / (1 + jnp.square(bin_centers) / jnp.square(d0))

    probs = jax.nn.softmax(logits, axis=-1)
    predicted_tm_term = (probs * tm_per_bin).sum(-1)

    cmap, pair_mask = _interchain_contact_weights(inputs, outputs, contact_dist)
    predicted_tm_term = predicted_tm_term * pair_mask

    normed = cmap / (1e-8 + cmap.sum(-1, keepdims=True))
    per_alignment = (predicted_tm_term * normed).sum(-1)
    return (per_alignment * residue_weights).max()


def add_actifptm_loss(self, weight=0.1, contact_dist=8.0):
    """actifpTM as a loss. NOTE: the .max() reduction means gradient reaches only one
    alignment frame -- the same sparsity BECraft criticises in ipTM. Good as a metric
    (weight=0); for a dense objective prefer add_actif_i_ptm_energy_loss below."""
    def loss_fn(inputs, outputs):
        return {"actifptm": 1 - get_actifptm(inputs, outputs, contact_dist)}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["actifptm"] = weight


def get_actif_iptmenergy(inputs, outputs, tm_lambda=1.0, contact_dist=8.0):
    """Hybrid: pTMEnergy's dense unnormalised energy, but weighted by contact
    probability like actifpTM instead of a uniform inter-chain mask. Keeps the dense
    gradient (no softmax, no max) while ignoring unbound flanks."""
    logits = outputs["predicted_aligned_error"]["logits"]
    breaks = outputs["predicted_aligned_error"]["breaks"]

    step = breaks[1] - breaks[0]
    bin_centers = jnp.append(breaks + step / 2, breaks[-1] + step / 2 + step)

    num_res = jnp.maximum(logits.shape[0], 19)
    d0 = 1.24 * (num_res - 15) ** (1.0 / 3) - 1.8
    g = 1.0 / (1 + jnp.square(bin_centers) / jnp.square(d0))

    weighted_logits = logits + tm_lambda * jnp.log(g)[None, None, :]
    positional_energy = -jax.nn.logsumexp(weighted_logits, axis=-1)

    w, _ = _interchain_contact_weights(inputs, outputs, contact_dist)
    return (positional_energy * w).sum() / jnp.maximum(w.sum(), 1.0)


def add_actif_i_ptm_energy_loss(self, weight=1.0, tm_lambda=1.0, contact_dist=8.0):
    def loss_fn(inputs, outputs):
        return {"actif_i_ptm_energy": get_actif_iptmenergy(inputs, outputs, tm_lambda, contact_dist)}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["actif_i_ptm_energy"] = weight

# add helicity loss
def add_helix_loss(self, weight=0):
    def binder_helicity(inputs, outputs):  
      if "offset" in inputs:
        offset = inputs["offset"]
      else:
        idx = inputs["residue_index"].flatten()
        offset = idx[:,None] - idx[None,:]

      # define distogram
      dgram = outputs["distogram"]["logits"]
      dgram_bins = get_dgram_bins(outputs)
      mask_2d = np.outer(np.append(np.zeros(self._target_len), np.ones(self._binder_len)), np.append(np.zeros(self._target_len), np.ones(self._binder_len)))

      x = _get_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
      if offset is None:
        if mask_2d is None:
          helix_loss = jnp.diagonal(x,3).mean()
        else:
          helix_loss = jnp.diagonal(x * mask_2d,3).sum() + (jnp.diagonal(mask_2d,3).sum() + 1e-8)
      else:
        mask = offset == 3
        if mask_2d is not None:
          mask = jnp.where(mask_2d,mask,0)
        helix_loss = jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

      return {"helix":helix_loss}
    self._callbacks["model"]["loss"].append(binder_helicity)
    self.opt["weights"]["helix"] = weight

# add N- and C-terminus distance loss
def add_termini_distance_loss(self, weight=0.1, threshold_distance=7.0):
    '''Add loss penalizing the distance between N and C termini'''
    def loss_fn(inputs, outputs):
        xyz = outputs["structure_module"]
        ca = xyz["final_atom_positions"][:, residue_constants.atom_order["CA"]]
        ca = ca[-self._binder_len:]  # Considering only the last _binder_len residues

        # Extract N-terminus (first CA atom) and C-terminus (last CA atom)
        n_terminus = ca[0]
        c_terminus = ca[-1]

        # Compute the distance between N and C termini
        termini_distance = jnp.linalg.norm(n_terminus - c_terminus)

        # Compute the deviation from the threshold distance using ELU activation
        deviation = jax.nn.elu(termini_distance - threshold_distance)

        # Ensure the loss is never lower than 0
        termini_distance_loss = jax.nn.relu(deviation)
        return {"NC": termini_distance_loss}

    # Append the loss function to the model callbacks
    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["NC"] = weight

# plot design trajectory losses
def plot_trajectory(af_model, design_name, design_paths):
    metrics_to_plot = ['loss', 'plddt', 'ptm', 'i_ptm', 'i_ptm_energy', 'actifptm', 'actif_i_ptm_energy', 'con', 'i_con', 'pae', 'i_pae', 'rg', 'mpnn']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    for index, metric in enumerate(metrics_to_plot):
        if metric in af_model.aux["log"]:
            # Create a new figure for each metric
            plt.figure()

            loss = af_model.get_loss(metric)
            # Create an x axis for iterations
            iterations = range(1, len(loss) + 1)

            plt.plot(iterations, loss, label=f'{metric}', color=colors[index % len(colors)])

            # Add labels and a legend
            plt.xlabel('Iterations')
            plt.ylabel(metric)
            plt.title(design_name)
            plt.legend()
            plt.grid(True)

            # Save the plot
            plt.savefig(os.path.join(design_paths["Trajectory/Plots"], design_name+"_"+metric+".png"), dpi=150)
            
            # Close the figure
            plt.close()