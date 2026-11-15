from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Any
import tempfile
import subprocess
import os
import re

import numpy as np

from Bio.PDB.Structure import Structure
from ..utils.structure_utils import get_single_chain,get_residues_by_chain
from ..utils.IO_utils import write_pdb
from ..utils.sequence_utils import normalize_aa_name_to_one_letter
from ..utils.logging_utils import Logger



def compute_pockets(struct: Structure, logger, min_rad: float = 1.8, max_rad: float =6.2, min_volume: int =50) -> List[Dict[str, Any]] | None:
    if min_rad <= 0 or max_rad <= 0 or min_volume <= 0 or min_rad > max_rad:
        logger.print("[ERROR] Invalid PyVOL parameters.")
        return None

    # ---------------------------
    # 1. 获取链和残基信息
    # ---------------------------
    chain = get_single_chain(struct, logger)
    if chain is None:
        return None

    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        return None

    # 建立 resseq -> resname 映射
    resseq_to_name: Dict[int, str] = {}
    for (hetflag, resseq, icode), resname, _ in residue_list:
        resseq_to_name[int(resseq)] = resname.upper()

    # ---------------------------
    # 2. 临时目录
    # ---------------------------
    with tempfile.TemporaryDirectory(prefix="pyvol_") as tmpdir:
        tmpdir = Path(tmpdir)

        pdb_path = tmpdir / "input.pdb"
        cfg_path = tmpdir / "pyvol.cfg"
        log_path = tmpdir / "pyvol.log"

        # ---------------------------
        # 3. 写 PDB
        # ---------------------------
        try:
            write_pdb(struct, pdb_path)
        except Exception as e:
            logger.print(f"[ERROR] Failed to write temporary PDB file for PyVOL: {e}")
            return None

        if not pdb_path.exists():
            logger.print(f"[ERROR] Temporary PDB file was not created.")
            return None

        # ---------------------------
        # 4. 写 PyVOL 配置
        # ---------------------------
        cfg_text = f"""[General]
prot_file = {pdb_path}
min_rad = {min_rad}
max_rad = {max_rad}
constrain_radii = True

[Specification]
mode = all
min_volume = {min_volume}

[Partitioning]
subdivide = False

[Output]
project_dir = {tmpdir}
prefix = run
logger_stream_level = INFO
logger_file_level = DEBUG
"""
        try:
            cfg_path.write_text(cfg_text, encoding="utf-8")
        except Exception as e:
            logger.print(f"[ERROR] Failed to write PyVOL config file: {e}")
            return None

        if not cfg_path.exists():
            logger.print(f"[ERROR] PyVOL config file was not created.")
            return None

        # ---------------------------
        # 5. 运行 PyVOL
        # ---------------------------
        env = os.environ.copy()
        env["TMPDIR"] = str(tmpdir)

        try:
            p = subprocess.run(
                ["pyvol", str(cfg_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                check=False,
            )
        except Exception as e:
            logger.print(f"[ERROR] Failed to run PyVOL: {e}")
            return None

        try:
            log_path.write_text(p.stdout or "", encoding="utf-8")
        except Exception as e:
            logger.print(f"[ERROR] Failed to write PyVOL log file: {e}")
            return None

        if p.returncode != 0:
            pyvol_output = (p.stdout or "").strip()
            if pyvol_output:
                pyvol_output = "\n".join(pyvol_output.splitlines()[-20:])
                logger.print(f"[ERROR] PyVOL failed with return code {p.returncode}. Output tail: {pyvol_output}")
            else:
                logger.print(f"[ERROR] PyVOL failed with return code {p.returncode}. This may be due to excessively small/large probe radius.")
            return None

        # ---------------------------
        # 6. 找 pocket 文件
        # ---------------------------
        obj_files = list(tmpdir.rglob("*_p*.obj"))
        xyzrg_files = list(tmpdir.rglob("*_p*.xyzrg"))

        rept_files = list(tmpdir.rglob("*.rept"))
        if not obj_files or not xyzrg_files:
            if rept_files:
                logger.print(f"[WARNING] No pockets detected by PyVOL.")
                return []
            else:
                logger.print(f"[ERROR] PyVOL output incomplete.")
                return None

        pocket_files: Dict[int, Dict[str, Path]] = {}
        pattern = re.compile(r"_p(\d+)\.", re.IGNORECASE)

        for pth in obj_files:
            m = pattern.search(pth.name)
            if m:
                pid = int(m.group(1))
                pocket_files.setdefault(pid, {})["obj"] = pth

        for pth in xyzrg_files:
            m = pattern.search(pth.name)
            if m:
                pid = int(m.group(1))
                pocket_files.setdefault(pid, {})["xyzrg"] = pth

        pocket_files = {
            pid: d for pid, d in pocket_files.items()
            if "obj" in d and "xyzrg" in d
        }

        if not pocket_files:
            logger.print(f"[ERROR] No valid pocket obj/xyzrg pairs found.")
            return None

        # ---------------------------
        # 7. 解析 volume（.rept）
        # ---------------------------
        volumes: Dict[int, float] = {}

        rept_files = list(tmpdir.rglob("*.rept"))
        if not rept_files:
            logger.print(f"[ERROR] No PyVOL report (.rept) files found.")
            return None

        csv_pid = re.compile(r"(?:_?p)(\d+)", re.IGNORECASE)

        for rp in rept_files:
            try:
                lines = rp.read_text(encoding="utf-8", errors="ignore").splitlines()
            except Exception as e:
                logger.print(f"[ERROR] Failed to read report file {rp}: {e}")
                return None

            if not lines:
                continue

            if lines[0].lower().startswith("name,volume"):
                for ln in lines[1:]:
                    parts = [x.strip() for x in ln.split(",")]
                    if len(parts) < 2:
                        continue
                    m = csv_pid.search(parts[0])
                    if not m:
                        continue
                    try:
                        volumes[int(m.group(1))] = float(parts[1])
                    except Exception:
                        continue

        # ---------------------------
        # 8. xyzrg -> 球 + residues + bbox
        # ---------------------------
        pocket_results: List[Dict[str, Any]] = []

        for pid, files in sorted(pocket_files.items()):
            xyzrg = files["xyzrg"]

            centers = []
            radii = []

            try:
                with open(xyzrg, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) < 4:
                            continue
                        try:
                            x, y, z, r = map(float, parts[:4])
                        except Exception:
                            continue
                        centers.append([x, y, z])
                        radii.append(r)
            except Exception as e:
                logger.print(f"[ERROR] Failed to read xyzrg file {xyzrg}: {e}")
                return None

            if not centers:
                logger.print(f"[ERROR] No valid spheres found in xyzrg file: {xyzrg}")
                return None

            centers = np.asarray(centers, dtype=float)
            radii = np.asarray(radii, dtype=float)

            # ---------------------------
            # 9. pocket center + box boundaries
            # 与输入 Structure 坐标系一致，因为 xyzrg 本身就是基于该结构坐标写出的
            # ---------------------------
            r = radii.reshape(-1, 1)
            mn = np.min(centers - r, axis=0)
            mx = np.max(centers + r, axis=0)

            pocket_center_coord = ((mn + mx) / 2.0).tolist()
            pocket_box_boundaries = (mx - mn).tolist()

            # ---------------------------
            # 10. 球 -> residues（最近CA）
            # ---------------------------
            residues_set = set()

            for center in centers:
                min_dist = float("inf")
                closest_res = None

                for (hetflag, resseq, icode), resname, coord in residue_list:
                    dx = coord[0] - center[0]
                    dy = coord[1] - center[1]
                    dz = coord[2] - center[2]
                    d2 = dx * dx + dy * dy + dz * dz

                    if d2 < min_dist:
                        min_dist = d2
                        closest_res = int(resseq)

                if closest_res is not None:
                    residues_set.add(closest_res)

            residues = []
            for resseq in sorted(residues_set):
                if resseq not in resseq_to_name:
                    continue
                residues.append(
                    {
                        "aa_id": int(resseq),
                        "aa_name": normalize_aa_name_to_one_letter(resseq_to_name[resseq]),
                    }
                )

            volume = volumes.get(pid)

            # 必须全部满足，否则跳过
            if (
                    volume is None
                    or len(centers) == 0
                    or not residues
                    or len(pocket_center_coord) != 3
                    or len(pocket_box_boundaries) != 3
            ):
                continue

            pocket_results.append(
                {
                    "volume": float(volume),
                    "n_spheres": int(len(centers)),
                    "residues": residues,
                    "pocket_center_coord": [float(x) for x in pocket_center_coord],
                    "pocket_box_boundaries": [float(x) for x in pocket_box_boundaries],
                }
            )

        if not pocket_results:
            logger.print(f"[WARNING] No valid pockets found after filtering.")
            return []

        # ---------------------------
        # 11. 按 volume 排序
        # ---------------------------
        pocket_results.sort(
            key=lambda x: (
                x["volume"] is not None,
                x["volume"] if x["volume"] is not None else -1.0,
            ),
            reverse=True,
        )

        return pocket_results

def calculate_pocket_statistics(pocket_regions: List[Dict[str, Any]]) -> Dict[str, Any]:

    pocket_num = len(pocket_regions)
    max_pocket_volume = 0.0
    total_pocket_volume = 0.0

    if pocket_num > 0:
        max_pocket_volume = max(float(pocket["volume"]) for pocket in pocket_regions)
        total_pocket_volume = sum(float(pocket["volume"]) for pocket in pocket_regions)

    return {
        "pocket_num": pocket_num,
        "max_pocket_volume": max_pocket_volume,
        "total_pocket_volume": total_pocket_volume,
    }



def postprocess_pocket_report_to_schema(
    raw_report: Dict[str, Any],
    logger: Logger | None = None,
) -> Dict[str, Any] | None:
    """
    Postprocess the raw EnzyWizard-Pocket report into the official JSON Schema field names.

    """

    if not isinstance(raw_report, dict):
        return None

    try:
        raw_statistics = raw_report.get("pocket_region_statistics")
        if not isinstance(raw_statistics, dict):
            return None

        raw_pocket_regions = raw_report.get("pocket_regions")
        if not isinstance(raw_pocket_regions, list):
            return None

        binding_pockets: List[Dict[str, Any]] = []

        for raw_pocket in raw_pocket_regions:
            if not isinstance(raw_pocket, dict):
                return None

            raw_residues = raw_pocket.get("residues")
            if not isinstance(raw_residues, list):
                return None

            residues: List[Dict[str, Any]] = []

            for raw_residue in raw_residues:
                if not isinstance(raw_residue, dict):
                    return None

                residues.append(
                    {
                        "residue_index": raw_residue.get("aa_id"),
                        "residue_name": raw_residue.get("aa_name"),
                    }
                )

            binding_pockets.append(
                {
                    "binding_pocket_volume": raw_pocket.get("volume"),
                    "binding_pocket_sphere_count": raw_pocket.get("n_spheres"),
                    "residues": residues,
                    "binding_pocket_center_coordinate": raw_pocket.get("pocket_center_coord"),
                    "binding_pocket_box_size": raw_pocket.get("pocket_box_boundaries"),
                }
            )

        schema_report: Dict[str, Any] = {
            "report_type": raw_report.get("output_type"),
            "binding_pocket_statistics": {
                "binding_pocket_count": raw_statistics.get("pocket_num"),
                "max_binding_pocket_volume": raw_statistics.get("max_pocket_volume"),
                "total_binding_pocket_volume": raw_statistics.get("total_pocket_volume"),
            },
            "binding_pockets": binding_pockets,
        }

        return schema_report

    except Exception as e:
        if logger is not None:
            logger.print(f"[ERROR] Failed to postprocess pocket report: {e}")
        return None

def generate_pocket_report(
    pocket_regions: List[Dict[str, Any]],
    logger: Logger | None = None,
) -> Dict[str, Any] | None:
    pocket_region_statistics = calculate_pocket_statistics(pocket_regions)

    raw_report = {
        "output_type": "enzywizard_pocket",
        "pocket_region_statistics": pocket_region_statistics,
        "pocket_regions": pocket_regions,
    }

    return postprocess_pocket_report_to_schema(raw_report, logger)
