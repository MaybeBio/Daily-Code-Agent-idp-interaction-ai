"""
3D visualization of the two-sphere phi(rho, z) profile produced by
spec/fisher_two_spheres.jl.

The spec saves the full 2D axisymmetric phi field per converged d under
res["profile_phi"][k] with companion res["profile_rho"][k],
res["profile_z"][k], res["profile_d"][k].  This script:

  1. Reads one (alpha, beta, d-index) tuple from the JSON output.
  2. Rotates the 2D phi(rho, z) field around the z-axis to recover the full
     3D axisymmetric volume.
  3. Renders the two spheres as semi-transparent meshes plus a translucent
     point cloud sampled with density proportional to (phi - phi_inf).  The
     visual effect is "denser cloud where phi is higher" -- wetting layers
     around each sphere show up as bright shells in the otherwise dilute box.

Usage:
    python viz_util/visualize_two_spheres_3d.py path/to/output_cpuK.json \\
        --alpha 1 --beta 1 --d_index 0 [--out figure.png]

NaN values in the JSON (JSON.jl encodes NaN as null -> Python None) are
converted in-place to np.nan when reading scalar arrays.
"""

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np


def _to_float(x):
    """Convert None -> NaN, otherwise float."""
    return float("nan") if x is None else float(x)


def _to_array(seq):
    return np.array([_to_float(v) for v in seq], dtype=float)


def load_result(path, alpha_key, beta_key):
    with open(path) as f:
        data = json.load(f)
    return data[str(alpha_key)][str(beta_key)]


def grid_density_cloud(phi, rho_grid, z_grid, R, z_c1, z_c2, phi_inf, n_theta):
    """
    Rotate the 2D (rho, z) grid around the z-axis using n_theta angular slices.
    Every grid point is included; alpha values are proportional to phi - phi_inf.
    Returns (xyz, phi_values, alpha_values).
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    # shape: (Nrho, Nz, n_theta)
    RHO  = rho_grid[:, np.newaxis, np.newaxis]
    ZZ   = z_grid[np.newaxis, :, np.newaxis]
    TH   = theta[np.newaxis, np.newaxis, :]
    PHI3 = phi[:, :, np.newaxis] * np.ones(n_theta)

    # x, y depend on (rho, theta) but must be broadcast across the z axis too,
    # so every array below is the full (Nrho, Nz, n_theta) shape.
    x = np.broadcast_to(RHO * np.cos(TH), PHI3.shape).ravel()
    y = np.broadcast_to(RHO * np.sin(TH), PHI3.shape).ravel()
    z = np.broadcast_to(ZZ, PHI3.shape).ravel()
    rho_flat = np.broadcast_to(RHO, PHI3.shape).ravel()
    phi_vals = PHI3.ravel()

    outside = ((rho_flat**2 + (z - z_c1)**2) > R**2) & \
              ((rho_flat**2 + (z - z_c2)**2) > R**2)
    x, y, z, phi_vals = x[outside], y[outside], z[outside], phi_vals[outside]

    phi_excess = np.clip(phi_vals - phi_inf, 0.0, None)
    amax = max(phi_excess.max(), 1e-12)
    # Drop the dilute bulk (alpha ~ 0): keep only the enriched wetting/bridge
    # shells so the scatter is small enough to render (the dense capillary
    # bridge of the new bridged branch shows up here as a bright neck).
    keep = phi_excess > 0.02 * amax
    x, y, z, phi_vals = x[keep], y[keep], z[keep], phi_vals[keep]
    alpha_vals = phi_excess[keep] / amax

    return np.column_stack([x, y, z]), phi_vals, alpha_vals


def render(result, d_index, out_path, n_theta):
    R       = float(result["R"])
    phi_inf = float(result["phi_inf"])
    omega_b = float(result["omega_b"])

    profile_d   = _to_array(result["profile_d"])
    if d_index < 0 or d_index >= len(profile_d):
        raise IndexError(f"d_index {d_index} out of range "
                         f"(have {len(profile_d)} converged d-points)")
    d   = float(profile_d[d_index])
    rho = _to_array(result["profile_rho"][d_index])
    z   = _to_array(result["profile_z"][d_index])
    phi = np.array(result["profile_phi"][d_index], dtype=float)  # (Nrho, Nz)

    # Label which branch this equilibrium profile came from (two-branch spec).
    d_grid_all = _to_array(result.get("d_grid", []))
    branch_arr = result.get("branch", None)
    branch_label = ""
    if branch_arr is not None and d_grid_all.size:
        kk = int(np.argmin(np.abs(d_grid_all - d)))
        branch_label = {1: ", branch=unbridged", 2: ", branch=bridged"}.get(
            int(branch_arr[kk]), "")

    z_c1 = 0.0
    z_c2 = 2.0 * R + d

    pts, phi_pts, alpha_pts = grid_density_cloud(
        phi, rho, z, R, z_c1, z_c2, phi_inf, n_theta=n_theta,
    )

    fig = plt.figure(figsize=(7.5, 11.0))
    ax  = fig.add_subplot(111, projection="3d")

    # Two spheres -- semi-opaque, light grey
    u  = np.linspace(0.0, 2.0 * np.pi, 48)
    v  = np.linspace(0.0, np.pi, 36)
    xs = R * np.outer(np.cos(u), np.sin(v))
    ys = R * np.outer(np.sin(u), np.sin(v))
    zs_c = np.outer(np.ones_like(u), np.cos(v))
    for zc in (z_c1, z_c2):
        ax.plot_surface(xs, ys, zc + R * zs_c,
                        color="#d0d0d0", alpha=0.55, linewidth=0,
                        edgecolor="none", shade=True, antialiased=True)

    # Density cloud -- alpha and color both encode phi
    norm = plt.Normalize(vmin=phi_inf, vmax=max(phi.max(), phi_inf + 1e-9))
    cmap = plt.cm.magma
    rgba = cmap(norm(phi_pts))
    rgba[:, 3] = alpha_pts
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
               c=rgba, s=5, edgecolor="none", depthshade=False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, shrink=0.55, pad=0.02)
    cb.set_label(r"$\phi$")

    # Equal aspect in (x, y); z is naturally longer than rho_max in this setup
    rho_max = rho[-1]
    z_span  = z[-1] - z[0]
    ax.set_box_aspect((1.0, 1.0, z_span / (2.0 * rho_max)))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(
        f"two-sphere wet field: phi_inf={phi_inf}, omega_b={omega_b}, "
        f"R={R}, d={d:.3f}{branch_label}"
        f"\nphi range = [{phi.min():.4f}, {phi.max():.4f}]"
    )
    ax.view_init(elev=0, azim=0)

    # Trim panes
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path",
                        help="Path to output_cpuK.json from spec/fisher_two_spheres.jl")
    parser.add_argument("--alpha", type=int, default=1,
                        help="Top-level alpha key (1-indexed)")
    parser.add_argument("--beta", type=int, default=1,
                        help="Inner beta key (1-indexed)")
    parser.add_argument("--d_index", type=int, default=0,
                        help="Index into res['profile_d'] (0-indexed; "
                             "only converged d-points are stored)")
    parser.add_argument("--n_theta", type=int, default=60,
                        help="Number of angular slices when rotating 2D grid to 3D")
    parser.add_argument("--out", default=None,
                        help="Output PNG path (default: alongside JSON)")
    args = parser.parse_args()

    result = load_result(args.json_path, args.alpha, args.beta)
    if args.out is None:
        base = os.path.splitext(os.path.basename(args.json_path))[0]
        outdir = os.path.dirname(os.path.abspath(args.json_path))
        args.out = os.path.join(outdir,
            f"{base}_3d_a{args.alpha}_b{args.beta}_d{args.d_index}.png")
    render(result, args.d_index, args.out, n_theta=args.n_theta)


if __name__ == "__main__":
    main()
