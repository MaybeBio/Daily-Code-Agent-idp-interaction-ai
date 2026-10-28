"""
2D cross-section visualization of the two-sphere phi(rho, z) profile produced by
spec/fisher_two_spheres.jl -- the rho-z companion to visualize_two_spheres_3d.py.

The spec saves, per converged separation d, the full 2D axisymmetric phi field in
res["profile_phi"][k] with res["profile_rho"][k], res["profile_z"][k],
res["profile_d"][k], plus res["branch"] (1=unbridged/dilute gap, 2=bridged/
capillary neck) aligned to res["d_grid"].  This script:

  1. Reads one (alpha, beta) result from the JSON output.
  2. Mirrors the axisymmetric phi(rho, z) field across the symmetry axis to a full
     -rho_max..+rho_max slice, recentred so the gap midplane sits at z=0 (both
     spheres straddle the centre).
  3. Renders it as a magma pcolormesh with the two spheres drawn as grey discs;
     cells inside either sphere are masked (the BVP only solves the exterior).

Each figure has two panels: the rho-z cross-section (left) and the gap-facing
surface concentration vs separation (right) -- phi on sphere 1's surface at the
pole pointing toward sphere 2 (rho=0, z=z_c1+R), plotted over all stored d with a
black dot marking the current separation.

Two modes:
  * single frame  -> PNG of one d-index (cross-section + curve with the dot at d).
  * --animate     -> GIF of the spheres approaching (capillary bridge forming)
                     then retreating (bridge breaking), with the dot sliding
                     along the surface-concentration curve.  Frames are sampled
                     n_outside from the unbridged branch + n_inside from the
                     bridged branch, ordered far->near, and played ping-pong
                     (forward then reversed) so one loop shows form-then-break.
                     The GIF and one vector PDF per distinct frame are written to
                     --outdir (created if missing).

Usage:
    # single cross-section
    python viz_util/visualize_two_spheres_2d.py path/to/output_cpuK.json \\
        --alpha 1 --beta 1 --d_index 0 [--out figure.png]

    # capillary-bridge animation (12 frames: 8 outside + 4 inside)
    python viz_util/visualize_two_spheres_2d.py path/to/output_cpuK.json \\
        --alpha 1 --beta 1 --animate --n_outside 8 --n_inside 4 \\
        --outdir media/figures/two_spheres_bridge

NaN values in the JSON (JSON.jl encodes NaN as null -> Python None) are converted
in-place to np.nan when reading scalar arrays.  The GIF is written with Pillow
(no imageio dependency); per-frame PDFs are vector output from matplotlib.
"""

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Circle
from PIL import Image


def _to_float(x):
    """Convert None -> NaN, otherwise float."""
    return float("nan") if x is None else float(x)


def _to_array(seq):
    return np.array([_to_float(v) for v in seq], dtype=float)


def load_result(path, alpha_key, beta_key):
    with open(path) as f:
        data = json.load(f)
    return data[str(alpha_key)][str(beta_key)]


def branch_of(result, d):
    """Branch label (1 unbridged, 2 bridged, 0 none) for separation d."""
    d_grid = _to_array(result.get("d_grid", []))
    branch = result.get("branch", None)
    if branch is None or d_grid.size == 0:
        return 0
    return int(branch[int(np.argmin(np.abs(d_grid - d)))])


def cross_section(result, d_index):
    """
    Build the mirrored, recentred rho-z cross-section for one stored d-index.

    Returns rho_signed (-rho_max..+rho_max), z_c (centred so the gap midplane is
    at 0), phi_full (masked NaN inside the spheres, shape (2*Nrho, Nz)), the two
    sphere centres (z, both on rho=0), R, and d.
    """
    R   = float(result["R"])
    d   = float(_to_array(result["profile_d"])[d_index])
    rho = _to_array(result["profile_rho"][d_index])
    z   = _to_array(result["profile_z"][d_index])
    phi = np.array(result["profile_phi"][d_index], dtype=float)  # (Nrho, Nz)

    # Spec geometry: sphere 1 at z=0, sphere 2 at z = 2R + d.  Recentre on the
    # gap midplane so the pair straddles z=0 symmetrically.
    z_mid  = R + 0.5 * d
    z_c    = z - z_mid
    c1, c2 = -(R + 0.5 * d), (R + 0.5 * d)   # sphere centres, symmetric about 0

    # Mask the interior of either sphere (junk values live there).
    RHO = rho[:, np.newaxis]
    Zc  = z_c[np.newaxis, :]
    inside = ((RHO**2 + (Zc - c1)**2) < R**2) | ((RHO**2 + (Zc - c2)**2) < R**2)
    phi = np.where(inside, np.nan, phi)

    # Mirror across the symmetry axis: rho -> -rho.
    rho_signed = np.concatenate([-rho[::-1], rho])
    phi_full   = np.concatenate([phi[::-1, :], phi], axis=0)  # (2*Nrho, Nz)

    return rho_signed, z_c, phi_full, (c1, c2), R, d


def draw_cross_section(ax, rho_signed, z_c, phi_full, centres, R, d, norm, cmap):
    """Render one cross-section onto ax; return the QuadMesh (for the colorbar)."""
    # Fill padding outside the (d-dependent) field extent with the dilute-bulk
    # colour so fixed-limit frames don't show white bands as the spheres move.
    ax.set_facecolor(cmap(0.0))
    mesh = ax.pcolormesh(rho_signed, z_c, phi_full.T,
                         cmap=cmap, norm=norm, shading="auto")
    for zc in centres:
        ax.add_patch(Circle((0.0, zc), R, facecolor="#d0d0d0",
                            edgecolor="#808080", linewidth=0.8, zorder=3))
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$z$")
    return mesh


def inner_pole_phi(result, d_index):
    """
    phi at the surface of sphere 1 on the pole that points toward sphere 2
    (along +z, into the gap): the node at rho=0, z = z_c1 + R = R (sphere 1 is
    centred at z=0 in the spec geometry).
    """
    R   = float(result["R"])
    rho = _to_array(result["profile_rho"][d_index])
    z   = _to_array(result["profile_z"][d_index])
    phi = np.array(result["profile_phi"][d_index], dtype=float)
    i0  = int(np.argmin(np.abs(rho - 0.0)))   # symmetry axis
    j   = int(np.argmin(np.abs(z - R)))       # gap-facing pole of sphere 1
    return float(phi[i0, j])


def surface_curve(result):
    """Gap-facing surface phi vs separation over all stored d (sorted by d)."""
    profile_d = _to_array(result["profile_d"])
    order = np.argsort(profile_d)
    d   = profile_d[order]
    phi = np.array([inner_pole_phi(result, int(k)) for k in order])
    br  = np.array([branch_of(result, dd) for dd in d])
    return d, phi, br


def draw_surface_curve(ax, d_curve, phi_curve, branch_curve, d_now, phi_now,
                       xlim=None, ylim=None):
    """phi(gap-facing pole) vs d, with a black dot at the current separation."""
    ax.plot(d_curve, phi_curve, "-", color="0.6", lw=1.3, zorder=1)
    for b, c, lab in ((1, "tab:blue", "unbridged"), (2, "tab:red", "bridged")):
        m = branch_curve == b
        if m.any():
            ax.plot(d_curve[m], phi_curve[m], "o", color=c, ms=5, zorder=2, label=lab)
    ax.axvline(d_now, color="0.85", lw=1.0, ls="--", zorder=0)
    ax.plot([d_now], [phi_now], "o", color="black", ms=11, zorder=5)
    ax.set_xlabel(r"separation $d$")
    ax.set_ylabel(r"surface $\phi$ at gap-facing pole")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(loc="center right", fontsize=8, framealpha=0.9)


def render_cross_section(result, d_index, out_path):
    rho_signed, z_c, phi_full, centres, R, d = cross_section(result, d_index)
    phi_inf = float(result["phi_inf"])
    omega_b = float(result["omega_b"])
    vmax    = max(np.nanmax(phi_full), phi_inf + 1e-9)
    norm    = plt.Normalize(vmin=phi_inf, vmax=vmax)

    d_curve, phi_curve, br_curve = surface_curve(result)
    phi_now = inner_pole_phi(result, d_index)
    lab = {1: "unbridged", 2: "bridged"}.get(branch_of(result, d), "?")

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(10.0, 8.5), gridspec_kw={"width_ratios": [1.0, 1.2]})
    mesh = draw_cross_section(axL, rho_signed, z_c, phi_full, centres, R, d,
                             norm, plt.cm.magma)
    cb = fig.colorbar(mesh, ax=axL, shrink=0.6, pad=0.02)
    cb.set_label(r"$\phi$")
    axL.set_title(f"cross-section: d={d:.3f}  ({lab})")
    draw_surface_curve(axR, d_curve, phi_curve, br_curve, d, phi_now)
    axR.set_title(f"gap-facing surface concentration\n"
                  f"phi_inf={phi_inf:.4f}, omega_b={omega_b}, R={R}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {out_path}")


def select_frames(result, n_outside, n_inside):
    """
    Pick stored d-indices: n_inside from the bridged branch + n_outside from the
    unbridged branch, evenly spaced, returned in approach order (descending d).
    """
    profile_d = _to_array(result["profile_d"])
    branches  = np.array([branch_of(result, d) for d in profile_d])
    inside  = np.where(branches == 2)[0]   # bridged   (small d)
    outside = np.where(branches == 1)[0]   # unbridged (large d)

    def pick(idx, n):
        if len(idx) == 0 or n <= 0:
            return np.array([], dtype=int)
        n = min(n, len(idx))
        return idx[np.unique(np.round(np.linspace(0, len(idx) - 1, n)).astype(int))]

    chosen = np.concatenate([pick(inside, n_inside), pick(outside, n_outside)])
    # Approach = far -> near = descending d.
    return chosen[np.argsort(-profile_d[chosen])]


def render_animation(result, gif_path, outdir, n_outside, n_inside, fps):
    frames_idx = select_frames(result, n_outside, n_inside)
    if len(frames_idx) < 2:
        raise ValueError("need >=2 frames; not enough stored branches in this result")

    secs = [cross_section(result, k) for k in frames_idx]

    phi_inf = float(result["phi_inf"])
    omega_b = float(result["omega_b"])
    vmax    = max(np.nanmax(s[2]) for s in secs)
    vmax    = max(vmax, phi_inf + 1e-9)
    norm    = plt.Normalize(vmin=phi_inf, vmax=vmax)
    cmap    = plt.cm.magma

    # Fixed cross-section axes so the spheres approach/retreat in a stable frame.
    max_d   = max(s[5] for s in secs)
    R       = secs[0][4]
    z_half  = (2.0 * R + 0.5 * max_d) * 1.05
    rho_max = secs[0][0][-1]

    # Fixed surface-curve axes (full stored d-range, with a little padding).
    d_curve, phi_curve, br_curve = surface_curve(result)
    dx = 0.04 * (d_curve.max() - d_curve.min())
    dy = max(phi_curve.max() - phi_curve.min(), 1e-3)
    cv_xlim = (d_curve.min() - dx, d_curve.max() + dx)
    cv_ylim = (phi_curve.min() - 0.18 * dy, phi_curve.max() + 0.18 * dy)

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(gif_path)) or ".", exist_ok=True)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
    })

    fig = plt.figure(figsize=(6.0, 5))
    canvas = FigureCanvasAgg(fig)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0])
    axL = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1])
    # Static colorbar (built once off a throwaway mappable).
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=axL, shrink=0.6, pad=0.02)
    cb.set_label(r"$\phi$")

    images = []
    for n, (rho_signed, z_c, phi_full, centres, R, d) in enumerate(secs):
        axL.clear(); axR.clear()
        draw_cross_section(axL, rho_signed, z_c, phi_full, centres, R, d, norm, cmap)
        axL.set_xlim(-rho_max, rho_max)
        axL.set_ylim(-z_half, z_half)
        lab = {1: "unbridged", 2: "bridged"}.get(branch_of(result, d), "?")
        axL.set_title(f"d={d:.3f}  ({lab})")
        phi_now = inner_pole_phi(result, int(frames_idx[n]))
        draw_surface_curve(axR, d_curve, phi_curve, br_curve, d, phi_now,
                           cv_xlim, cv_ylim)
        axR.set_title(f"gap-facing surface concentration\n"
                      f"phi_inf={phi_inf:.4f}, omega_b={omega_b}, R={R}")
        fig.tight_layout()
        # Per-frame vector PDF (the 12 distinct far->near frames).
        fig.savefig(os.path.join(outdir, f"frame_{n:02d}_d{d:05.2f}.pdf"),
                    bbox_inches="tight")
        # Raster snapshot (fixed canvas size) for the GIF.
        canvas.draw()
        images.append(Image.fromarray(np.asarray(canvas.buffer_rgba())).convert("RGB"))
    plt.close(fig)

    # Ping-pong: approach then retreat (drop both endpoints on the way back so no
    # frame is duplicated at the turnarounds), so one loop shows form-then-break.
    sequence = images + images[-2:0:-1]
    duration = int(round(1000.0 / max(fps, 1e-6)))
    sequence[0].save(gif_path, save_all=True, append_images=sequence[1:],
                     duration=duration, loop=0, disposal=2)
    ds = [round(float(s[5]), 2) for s in secs]
    print(f"saved -> {gif_path}  ({len(images)} frames far->near at d={ds}, "
          f"ping-pong loop)")
    print(f"saved -> {outdir}/frame_*.pdf  ({len(images)} vector frames)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("json_path",
                        help="Path to output_cpuK.json from spec/fisher_two_spheres.jl")
    parser.add_argument("--alpha", type=int, default=1, help="Top-level alpha key (1-indexed)")
    parser.add_argument("--beta", type=int, default=1, help="Inner beta key (1-indexed)")
    parser.add_argument("--d_index", type=int, default=0,
                        help="Index into res['profile_d'] (single-frame mode)")
    parser.add_argument("--animate", action="store_true",
                        help="Build the approach/retreat capillary-bridge GIF")
    parser.add_argument("--n_outside", type=int, default=8,
                        help="Frames sampled from the unbridged branch (large d)")
    parser.add_argument("--n_inside", type=int, default=4,
                        help="Frames sampled from the bridged branch (small d)")
    parser.add_argument("--fps", type=float, default=3.0, help="GIF frames per second")
    parser.add_argument("--outdir", default=None,
                        help="Directory for the animation's GIF + per-frame PDFs "
                             "(created if missing; default beside JSON)")
    parser.add_argument("--out", default=None,
                        help="Explicit output path (single-frame PNG, or GIF path "
                             "overriding the default inside --outdir)")
    args = parser.parse_args()

    result = load_result(args.json_path, args.alpha, args.beta)
    base   = os.path.splitext(os.path.basename(args.json_path))[0]
    jdir   = os.path.dirname(os.path.abspath(args.json_path))
    tag    = f"{base}_2d_bridge_a{args.alpha}_b{args.beta}"

    if args.animate:
        outdir   = args.outdir or os.path.join(jdir, tag)
        gif_path = args.out or os.path.join(outdir, f"{tag}.gif")
        render_animation(result, gif_path, outdir,
                         args.n_outside, args.n_inside, args.fps)
    else:
        out = args.out or os.path.join(
            jdir, f"{base}_2d_a{args.alpha}_b{args.beta}_d{args.d_index}.png")
        render_cross_section(result, args.d_index, out)


if __name__ == "__main__":
    main()
