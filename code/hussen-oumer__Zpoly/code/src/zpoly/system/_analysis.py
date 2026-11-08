"""Tg and density from the cooling steps of a bulk system"""

# rq-a5952119

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from ._errors import SystemAnalysisError

MIN_POINTS = 3


@dataclass(frozen=True)
class TgResult:
    tg: float
    density_300: float | None
    temperatures: tuple[float, ...]
    densities: tuple[float, ...]
    std: tuple[float, ...]
    fit: dict


# rq-a5952119
def fit_tg(temperatures, densities) -> dict:
    """Bilinear fit"""
    t = np.asarray(temperatures, dtype=float)
    d = np.asarray(densities, dtype=float)
    order = np.argsort(t)
    t, d = t[order], d[order]
    if len(t) < 2 * MIN_POINTS:
        raise SystemAnalysisError(f"A Tg fit needs at least {2 * MIN_POINTS} temperatures, got {len(t)}.")
    best = None
    for k in range(MIN_POINTS, len(t) - MIN_POINTS + 1):
        lo = np.polyfit(t[:k], d[:k], 1)
        hi = np.polyfit(t[k:], d[k:], 1)
        sse = float(((np.polyval(lo, t[:k]) - d[:k]) ** 2).sum() + ((np.polyval(hi, t[k:]) - d[k:]) ** 2).sum())
        if best is None or sse < best[0]:
            best = (sse, lo, hi)
    sse, lo, hi = best
    if abs(lo[0] - hi[0]) < 1e-15:
        raise SystemAnalysisError("The two fitted lines are parallel; no Tg can be determined.")
    tg = float((hi[1] - lo[1]) / (lo[0] - hi[0]))
    return {"tg": tg, "glass": {"slope": float(lo[0]), "intercept": float(lo[1])},
            "melt": {"slope": float(hi[0]), "intercept": float(hi[1])}, "sse": sse}


def _density(gmx: str, edr: Path, begin_ps: float) -> tuple[float, float]:
    xvg = edr.with_name(edr.stem + "_density.xvg")
    try:
        proc = subprocess.run([gmx, "energy", "-f", str(edr), "-o", str(xvg), "-b", f"{begin_ps:g}"],
                              input="Density\n\n", capture_output=True, text=True, timeout=600)
    except FileNotFoundError:
        raise SystemAnalysisError(f"'{gmx}' was not found; install GROMACS or pass gmx=...") from None
    if proc.returncode != 0 or not xvg.is_file():
        raise SystemAnalysisError(f"gmx energy failed for {edr.name}:\n" + (proc.stdout + proc.stderr)[-1500:])
    values = [float(line.split()[1]) for line in xvg.read_text().splitlines()
              if line.strip() and not line.startswith(("#", "@"))]
    if not values:
        raise SystemAnalysisError(f"No density data in {xvg.name}.")
    arr = np.array(values) / 1000.0  # kg/m3 -> g/cm3
    return float(arr.mean()), float(arr.std())


# rq-a5952119
def analyze_tg(folder: str | os.PathLike, *, average_last_ns: float = 1.0, gmx: str = "gmx") -> TgResult:
    base = Path(folder)
    try:
        protocol = yaml.safe_load((base / "protocol.yaml").read_text(encoding="utf-8"))
    except OSError as exc:
        raise SystemAnalysisError(f"{base} has no protocol.yaml (was it written by bulk_system?)") from exc
    cool = protocol["cooling"]
    exe = shutil.which(gmx) or gmx
    temps, dens, stds = [], [], []
    for step in cool["steps"]:
        edr = base / f"{step['name']}.edr"
        if not edr.is_file():
            raise SystemAnalysisError(f"{step['name']}: {edr.name} is missing (has the protocol finished?)")
        begin = max(0.0, (step["time_ns"] - average_last_ns) * 1000.0)
        mean, std = _density(exe, edr, begin)
        temps.append(float(step["temperature"]))
        dens.append(mean)
        stds.append(std)
    fit = fit_tg(temps, dens)
    if 300.0 in temps:
        d300 = dens[temps.index(300.0)]
    elif min(temps) < 300.0 < max(temps):
        order = np.argsort(temps)
        d300 = float(np.interp(300.0, np.array(temps)[order], np.array(dens)[order]))
    else:
        d300 = None
    rows = zip(temps, dens, stds, strict=True)
    lines = ["temperature_K,density_g_cm3,std_g_cm3"] + [f"{t:g},{d:.5f},{s:.5f}" for t, d, s in rows]
    (base / "tg_density.csv").write_text("\n".join(lines) + "\n", encoding="ascii")
    summary = {"tg_K": round(fit["tg"], 2), "density_300K_g_cm3": None if d300 is None else round(d300, 5),
               "fit": fit, "average_last_ns": average_last_ns}
    (base / "tg_result.yaml").write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
    d300_txt = "n/a" if d300 is None else f"{d300:.4f} g/cm3"
    print(f"Tg = {fit['tg']:.1f} K, density at 300 K = {d300_txt} ({len(temps)} temperatures) -> "
          f"{base / 'tg_density.csv'}, {base / 'tg_result.yaml'}")
    return TgResult(fit["tg"], d300, tuple(temps), tuple(dens), tuple(stds), fit)


# LLPS slab
DILUTE_GAP_NM = 1.5


@dataclass(frozen=True)
class SlabResult:
    temperature: float
    dense: float | None
    dilute: float | None
    width_nm: float
    profile: tuple


# rq-1dce521e
def phase_concentrations(z, rho) -> tuple[float, float | None, float]:
    """Dense, dilute and width from a density profile."""
    z = np.asarray(z, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if len(z) < 5:
        raise SystemAnalysisError("The density profile has too few bins.")
    dz = float(np.median(np.diff(z)))
    top = int(np.argmax(rho))
    half = 0.5 * rho[top]
    left = top
    while left > 0 and rho[left - 1] >= half:
        left -= 1
    right = top
    while right < len(rho) - 1 and rho[right + 1] >= half:
        right += 1
    width = (right - left + 1) * dz
    n = right - left + 1
    core = rho[left + n // 4: right - n // 4 + 1] if n >= 4 else rho[left:right + 1]
    dense = float(core.mean())
    z_left, z_right = z[left], z[right]
    far = [(r) for zz, r in zip(z, rho, strict=True)
           if (zz < z_left - DILUTE_GAP_NM) or (zz > z_right + DILUTE_GAP_NM)]
    dilute = float(np.mean(far)) if far else None
    return dense, dilute, float(width)


def _profile(gmx: str, base: Path, name: str, slices: int, skip_ps: float):
    xvg = base / f"{name}_density.xvg"
    try:
        proc = subprocess.run([gmx, "density", "-s", str(base / f"{name}.tpr"), "-f", str(base / f"{name}.xtc"),
                               "-n", str(base / "index.ndx"), "-o", str(xvg), "-dens", "mass", "-center",
                               "-sl", str(slices), "-b", f"{skip_ps:g}"],
                              input="Polymer\nPolymer\n", capture_output=True, text=True, timeout=1800)
    except FileNotFoundError:
        raise SystemAnalysisError(f"'{gmx}' was not found; install GROMACS or pass gmx=...") from None
    if proc.returncode != 0 or not xvg.is_file():
        raise SystemAnalysisError(f"gmx density failed for {name}:\n" + (proc.stdout + proc.stderr)[-1500:])
    rows = [line.split()[:2] for line in xvg.read_text().splitlines()
            if line.strip() and not line.startswith(("#", "@"))]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


# rq-1dce521e
def analyze_slab(folder: str | os.PathLike, *, skip_ns: float = 50.0, bin_nm: float = 0.1,
                 gmx: str = "gmx") -> list[SlabResult]:
    """Phase concentrations per production temperature."""
    import warnings

    base = Path(folder)
    try:
        protocol = yaml.safe_load((base / "protocol.yaml").read_text(encoding="utf-8"))
    except OSError as exc:
        raise SystemAnalysisError(f"{base} has no protocol.yaml (was it written by slab_system?)") from exc
    slices = max(5, round(protocol["box_nm"][2] / bin_nm))
    exe = shutil.which(gmx) or gmx
    results = []
    for prod in protocol["productions"]:
        name = prod["name"]
        if not (base / f"{name}.xtc").is_file():
            raise SystemAnalysisError(f"{name}: {name}.xtc is missing (has this temperature been run?)")
        z, rho = _profile(exe, base, name, slices, skip_ns * 1000.0)
        dense, dilute, width = phase_concentrations(z, rho)
        if dilute is None:
            warnings.warn(f"{name}: no dilute region more than {DILUTE_GAP_NM} nm from the slab; the box is too "
                          "short for a dilute-phase concentration.", stacklevel=2)
        (base / f"density_profile_{name.removeprefix('prod_')}.csv").write_text(
            "z_nm,polymer_density_g_L\n" + "".join(f"{a:.4f},{b:.4f}\n" for a, b in zip(z, rho, strict=True)),
            encoding="ascii")
        results.append(SlabResult(float(prod["temperature"]), dense, dilute, width, tuple(zip(z, rho, strict=True))))
    lines = ["temperature_K,dense_g_L,dilute_g_L,slab_width_nm"]
    for r in results:
        dil = "" if r.dilute is None else f"{r.dilute:.3f}"
        lines.append(f"{r.temperature:g},{r.dense:.3f},{dil},{r.width_nm:.3f}")
    (base / "phase_diagram.csv").write_text("\n".join(lines) + "\n", encoding="ascii")
    for r in results:
        dil = "n/a" if r.dilute is None else f"{r.dilute:.1f}"
        print(f"{r.temperature:g} K: dense {r.dense:.1f} g/L, dilute {dil} g/L, slab {r.width_nm:.2f} nm")
    print(f"-> {base / 'phase_diagram.csv'}")
    return results
