"""
Compute bound and unbound times using hysteresis definition.

Replaces the binary threshold approach (is_bound = dist < 4) with a
hysteresis-based definition that prevents chatter in the gap region.

Definition:
    t^0_o = 0
    t^i_x = min{t > t^i_o : d_t > d_o}
    t^{i+1}_o = min{t > t^i_x : d_t < d_x}
    T^i_o = t^i_x - t^i_o  (unbound times)
    T^i_x = t^{i+1}_o - t^i_x  (bound times)

With d_o = 8 (separation threshold) and d_x = 2 (contact threshold).
"""

import numpy as np


def compute_hysteresis_times(pattern_distances, d_o=8.0, d_x=2.0, nruns=16, nparams=16):
    """
    Compute bound and unbound times using hysteresis thresholds.

    Args:
        pattern_distances: 16×16 list of lists of numpy arrays (inhomogeneous length)
        d_o: outer/separation threshold (default 8.0) — system is unbound when d > d_o
        d_x: inner/contact threshold (default 2.0) — system is bound when d < d_x

    Returns:
        (bound_times, unbound_times): each is 16×16 list of lists of durations
            bound_times[i][j]   = list of T^i_x durations (time to unbind from contact)
            unbound_times[i][j] = list of T^i_o durations (time to rebind after separation)
    """
    bound_times = [[[] for _ in range(nruns)] for _ in range(nparams)]
    unbound_times = [[[] for _ in range(nruns)] for _ in range(nparams)]
    
    for i in range(nparams):
        for j in range(nruns):
            traj = np.asarray(pattern_distances[i][j])
            if len(traj) == 0:
                continue

            # Find all indices where d > d_o (separation points)
            sep_idx = np.where(traj > d_o)[0]
            # Find all indices where d < d_x (contact points)
            con_idx = np.where(traj < d_x)[0]

            if len(sep_idx) == 0 or len(con_idx) == 0:
                continue

            idx_o = 0
            for sep_i in sep_idx:
                if sep_i < idx_o:
                    continue

                # Record unbound duration
                T_o = sep_i - idx_o
                unbound_times[i][j].append(T_o)

                # Find next contact after this separation
                next_con = con_idx[con_idx > sep_i]
                if len(next_con) == 0:
                    break

                idx_o_next = next_con[0]
                T_x = idx_o_next - sep_i
                bound_times[i][j].append(T_x)

                idx_o = idx_o_next

    return bound_times, unbound_times
