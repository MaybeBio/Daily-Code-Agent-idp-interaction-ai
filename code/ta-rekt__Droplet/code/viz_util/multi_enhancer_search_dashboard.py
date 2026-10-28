def plot_simulation_dashboard(res_or_path, keya, keyb, show=True, psize=30, diff_window=1000):
    """
    Plot a comprehensive simulation dashboard with:
    - 0,0: Acceptance rate
    - 0,1: Max pairwise distance between patterns
    - 1,0: Hamiltonians (poly and tf)
    - 1,1: Effective diffusion coefficients for each pattern
    - 2,0: Interaction matrix heatmap (mainPatterns vs tf only)
    - 2,1: Max Droplet Size
    - 3,: Final configuration (3D scatter)
    
    res_or_path: dict or path to JSON file produced by the Julia simulation.
    diff_window: window size (frames) for computing diffusion coefficients
    """
    # Load result
    if isinstance(res_or_path, str):
        with open(res_or_path, 'r') as f:
            res_ = json.load(f)
            res = res_[str(keya)][str(keyb)]
    else:
        res = res_or_path[str(keya)][str(keyb)]

    # Safely extract arrays
    def to_np(key):
        arr = np.array(res.get(key, []))
        return arr
    
    intmat = to_np("interaction_matrix")
    acc = to_np("acceptance_rate")
    pattern_dists_raw = res.get("pattern_dists", [])  # list of dicts with pairwise distances
    pattern_centers_raw = res.get("pattern_centers", [])  # list of 3x5 matrices
    droplet_size = to_np("max_droplet_size")
    ham_poly1 = to_np("hamiltonian_poly")
    ham_poly2 = to_np("hamiltonian_tf")

    print("Species name: ", res.get("species_name", "N/A"))
    # tf self-interaction is at index corresponding to tf (last species)
    tf_idx = intmat.shape[0] - 1 if intmat.size > 0 else 4
    print("Temperature: ", -1.0 / round(intmat[tf_idx, tf_idx], 2) if intmat.size > 0 else "N/A")
    print("Volume fraction: ", res["length_tf"]/res["boxlength"]**3)
    
    # Process pattern distances: extract max pairwise distance at each frame
    max_pattern_dists = []
    for dist_dict in pattern_dists_raw:
        if isinstance(dist_dict, dict) and len(dist_dict) > 0:
            max_pattern_dists.append(max(dist_dict.values()))
        else:
            max_pattern_dists.append(np.nan)
    max_pattern_dists = np.array(max_pattern_dists)
    
    # Process pattern centers: convert to array of shape (n_frames, 3, 5)
    pattern_centers = []
    for centers in pattern_centers_raw:
        if isinstance(centers, (list, np.ndarray)):
            pattern_centers.append(np.array(centers))
    pattern_centers = np.array(pattern_centers) if len(pattern_centers) > 0 else np.array([])
    
    # Compute effective diffusion coefficients for each pattern using MSD
    def compute_diffusion_coeffs(centers, window):
        """
        Compute effective diffusion coefficient D from MSD over a sliding window.
        D = MSD / (6 * dt) for 3D diffusion, where dt = window size
        Returns array of shape (n_windows, 5) with D for each pattern
        """
        if len(centers) < window:
            return np.array([]), np.array([])
        
        n_frames = len(centers)
        n_windows = n_frames - window
        diff_coeffs = np.zeros((n_windows, 5))
        window_starts = np.arange(n_windows)
        
        for i in range(n_windows):
            start_pos = centers[i]  # shape (3, 5)
            end_pos = centers[i + window]  # shape (3, 5)
            
            # MSD for each pattern (sum over x, y, z)
            msd = np.sum((end_pos - start_pos)**2, axis=0)  # shape (5,)
            
            # D = MSD / (6 * dt), dt = window
            diff_coeffs[i] = msd / (6 * window)
        
        return window_starts, diff_coeffs
    
    # Determine number of frames from the longest array
    n = max(len(acc), len(max_pattern_dists), len(ham_poly1), len(ham_poly2), len(droplet_size))
    frames = np.arange(1, n + 1)

    mono_colors = ["gray", "yellow", "red", "yellow", 
                   "yellow", "yellow", "red", "yellow"]

    # Create figure with 4x2 gridspec
    fig = plt.figure(figsize=(8,14))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 2], hspace=0.35, wspace=0.3)

    # ========== Row 0: Acceptance rate and Max pairwise distance ==========
    
    # (0, 0) Acceptance rate
    ax_acc = fig.add_subplot(gs[0, 0])
    if len(acc) > 1:
        ax_acc.plot(frames[:len(acc)], acc, color='C0', linewidth=1.5)
    ax_acc.set_title('Acceptance Rate', fontsize=12, fontweight='bold')
    ax_acc.set_xlabel('Frame')
    ax_acc.set_ylabel('Rate')
    ax_acc.grid(True, alpha=0.3)

    # (0, 1) Max pairwise distance between patterns
    ax_max_dist = fig.add_subplot(gs[0, 1])
    if len(max_pattern_dists) > 1:
        ax_max_dist.plot(np.arange(1, len(max_pattern_dists) + 1), max_pattern_dists, 
                         color='C1', linewidth=1.5)
    ax_max_dist.set_title('Max Pairwise Pattern Distance', fontsize=12, fontweight='bold')
    ax_max_dist.set_xlabel('Frame')
    ax_max_dist.set_ylabel('Distance')
    ax_max_dist.grid(True, alpha=0.3)

    # ========== Row 1: Hamiltonians and Diffusion coefficients ==========
    
    # (1, 0) Hamiltonians (poly and tf on same plot)
    ax_ham_poly = fig.add_subplot(gs[1, 0])
    
    color_poly1 = 'C5'
    color_poly2 = 'C9'
    ax_ham_poly.set_xlabel('Frame')
    ax_ham_poly.set_ylabel('Polymer Hamiltonian', color=color_poly1)
    if len(ham_poly1) > 1:
        ax_ham_poly.plot(frames[:len(ham_poly1)], ham_poly1, color=color_poly1, 
                        linewidth=1.5, label='Poly')
    ax_ham_poly.tick_params(axis='y', labelcolor=color_poly1)
    ax_ham_poly.grid(True, alpha=0.3)
    
    # Create second y-axis for tf hamiltonian
    ax_ham_tf = ax_ham_poly.twinx()
    color_tf = 'C6'
    ax_ham_tf.set_ylabel('TF Hamiltonian', color=color_tf)
    if len(ham_poly2) > 1:
        ax_ham_tf.plot(frames[:len(ham_poly2)], ham_poly2, color=color_tf, 
                        linewidth=1.5, linestyle='--', label='TF')
    ax_ham_tf.tick_params(axis='y', labelcolor=color_tf)
    
    # Add title and legend
    ax_ham_poly.set_title('Hamiltonians', fontsize=12, fontweight='bold')
    lines1, labels1 = ax_ham_poly.get_legend_handles_labels()
    lines2, labels2 = ax_ham_tf.get_legend_handles_labels()
    ax_ham_poly.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

    # (1, 1) Effective diffusion coefficients for each pattern
    ax_diff = fig.add_subplot(gs[1, 1])
    if len(pattern_centers) > 0:
        window_starts, diff_coeffs = compute_diffusion_coeffs(pattern_centers, diff_window)
        if len(diff_coeffs) > 0:
            pattern_labels = ['P1', 'P2', 'P3', 'P4', 'P5']
            colors = ['C0', 'C1', 'C2', 'C3', 'C4']
            for i in range(5):
                ax_diff.plot(window_starts, diff_coeffs[:, i], color=colors[i], 
                            linewidth=1.5, label=pattern_labels[i], alpha=0.8)
            ax_diff.legend(fontsize=8, ncol=2)
    ax_diff.set_title(f'Effective Diffusion (window={diff_window})', fontsize=12, fontweight='bold')
    ax_diff.set_xlabel('Frame')
    ax_diff.set_ylabel('D (a.u.)')
    ax_diff.grid(True, alpha=0.3)

    # ========== Row 2: Interaction matrix and Max Droplet Size ==========
    
    # (2, 0) Interaction matrix heatmap (mainPatterns vs tf only)
    ax_intmat = fig.add_subplot(gs[2, 0])
    
    # Extract submatrix: mainPattern indices are 3,5,7,9,11 (0-indexed: 2,4,6,8,10)
    # tf index is the last one (index 12 -> 0-indexed: 11 or shape-1)
    if intmat.size > 0:
        # mainPattern indices (Julia 1-indexed: 3,5,7,9,11 -> Python 0-indexed: 2,4,6,8,10)
        pattern_indices = [2, 4, 6, 8, 10]
        tf_idx = intmat.shape[0] - 1  # last index is tf
        
        # Create submatrix of pattern-tf interactions
        sub_indices = pattern_indices + [tf_idx]
        submat = intmat[np.ix_(sub_indices, sub_indices)]
        
        labels = ["P1", "P2", "P3", "P4", "P5", "TF"]
        
        maxv = np.max(np.abs(submat)) if np.max(np.abs(submat)) > 0 else 1
        norm = mpl.colors.Normalize(vmin=-maxv, vmax=maxv)
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")
        im = ax_intmat.matshow(submat, norm=norm, cmap="coolwarm", alpha=0.7)
        ax_intmat.set_xticks(range(len(labels)))
        ax_intmat.set_xticklabels(labels, rotation=45, ha="left")
        ax_intmat.set_yticks(range(len(labels)))
        ax_intmat.set_yticklabels(labels, rotation=45, va="top")
        ax_intmat.set_title('Pattern-TF Interactions', fontsize=12, fontweight='bold', pad=20)
        cbar = plt.colorbar(mappable, ax=ax_intmat, shrink=0.8)
        cbar.set_label('Interaction energy')

    # (2, 1) Max Droplet Size
    ax_size = fig.add_subplot(gs[2, 1])
    if len(droplet_size) > 1:
        ax_size.plot(frames[:len(droplet_size)], droplet_size, color='C2', linewidth=1.5)
    ax_size.set_title('Max Droplet Size', fontsize=12, fontweight='bold')
    ax_size.set_xlabel('Frame')
    ax_size.set_ylabel('Size (beads)')
    ax_size.grid(True, alpha=0.3)

    # ========== Row 3: Final configuration (3D scatter) ==========
    
    # 3D final configuration
    ax3d = fig.add_subplot(gs[3, :], projection='3d')
    
    # Extract polymer coordinates (transposed from Julia: cols are x, y, z, type)
    poly_coords = np.array(res.get("final_coordinates_poly", []))
    tf_coords = np.array(res.get("final_coordinates_tf", []))
    
    if poly_coords.size > 0:
        # Handle both row-major and column-major layouts
        if poly_coords.shape[0] == 4:  # shape is (4, N)
            x1, y1, z1, type1 = poly_coords
        else:  # shape is (N, 4)
            x1, y1, z1, type1 = poly_coords.T
        x1 = np.array(x1); y1 = np.array(y1); z1 = np.array(z1); type1 = np.array(type1, dtype=int)
        n1 = min(len(x1), len(type1))
        ax3d.scatter(x1[:n1], y1[:n1], z1[:n1], s=psize,
                     c=np.array(mono_colors)[type1[:n1]-1], marker='o',
                     edgecolor='none', alpha=1)
    
    if tf_coords.size > 0:
        # Handle both row-major and column-major layouts
        if tf_coords.shape[0] == 4:  # shape is (4, N)
            xb, yb, zb, _ = tf_coords
        else:  # shape is (N, 4)
            xb, yb, zb, _ = tf_coords.T
        xb = np.array(xb); yb = np.array(yb); zb = np.array(zb)
        ax3d.scatter(xb, yb, zb, s=psize, c='blue', marker='o',
                     edgecolor='none', alpha=0.1)

    L = res.get('boxlength', None)
    if L is not None:
        ax3d.set_xlim(0, L)
        ax3d.set_ylim(0, L)
        ax3d.set_zlim(0, L)
    ax3d.set_title('Final polymer configuration')
    ax3d.set_xlabel('x')
    ax3d.set_ylabel('y')
    ax3d.set_zlabel('z')

    # Use subplots_adjust instead of tight_layout to avoid warnings with 3D and twinx axes
    plt.subplots_adjust(hspace=0.35, wspace=0.3, top=0.96, bottom=0.05, left=0.08, right=0.95)
    if show:
        plt.show()
    return fig, res