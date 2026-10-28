import numpy as np
import mdtraj as md
import scipy

def kd_calculation(pdb,number_contact,weights):
    w_eq = np.array(weights, dtype=float)

    traj = md.load(pdb)
    box_lengths = traj.unitcell_lengths  

    box_lengths = traj.unitcell_lengths  
    v = np.prod(box_lengths, axis=1)  # volume per frame in nm^3
    volume_L = v.mean() * 1e-24       # average volume in liters

    # Load and normalize weights
    # w_file = "/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/reweight/G5/gabis_weights/weights_corr_G5"
    # w_eq = np.loadtxt(w_file)
    # w_eq /= w_eq.sum()

    bound_mask = np.array(number_contact) != 0
    unbound_mask = ~bound_mask

    w_bound = w_eq[bound_mask].sum()
    w_unbound = w_eq[unbound_mask].sum()

    Kd_pop_weighted = (w_unbound / w_bound) * (1 / (scipy.constants.N_A * volume_L))
    #print(f"Kd from populations: {Kd_pop_weighted:.3e} M")

    frame_time_ns = 0.02
    frame_time_s = frame_time_ns * 1e-9

    T_bound_s = (w_eq[bound_mask].sum()) * frame_time_s
    T_unbound_s = (w_eq[unbound_mask].sum()) * frame_time_s

    # Count transitions (these stay unweighted)
    n_binding = sum((number_contact[i-1] == 0 and number_contact[i] != 0) for i in range(1, len(number_contact)))
    n_unbinding = sum((number_contact[i-1] != 0 and number_contact[i] == 0) for i in range(1, len(number_contact)))

    # Use known ligand concentration (or compute from box volume)
    ligand_concentration = 1 / (scipy.constants.N_A * volume_L)  # mol/L

    # Kinetic rates
    k_off = n_unbinding / T_bound_s
    k_on = n_binding / (T_unbound_s * ligand_concentration)

    Kd_kinetic_weighted = k_off / k_on
    #print(f"Kd from kinetics: {Kd_kinetic_weighted:.3e} M")

    return Kd_kinetic_weighted, Kd_pop_weighted


