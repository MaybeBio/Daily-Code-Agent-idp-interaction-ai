#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_global_contact_log.py

- 读取 extracted_500ns_bulk_chains 下的轨迹文件: chainXXX_bulk_500ns.dcd
- 计算系综平均 Contact Map
- LogNorm 绘图
- 输出 Sum of P_ij (上三角所有 pair 的概率和)
"""

import os
import re
import glob
import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from tqdm import tqdm

# ==================== 配置区域 ====================
TRAJ_DIR = "./extracted_500ns_bulk_chains"
PDB_FILE = "starts/seed_000.pdb"

RCUT = 0.66
GAP  = 2

VMIN = 1e-4
VMAX = 1e-2
STRIDE = 10
OUT_DIR = "analysis_results_bulk"
CHUNK = 1000  # iterload chunk size
# =================================================

def natural_key_chain(path):
    """
    chain000_bulk_500ns.dcd -> 0
    用于按 chain 编号排序
    """
    base = os.path.basename(path)
    m = re.search(r"chain(\d+)_bulk_500ns\.dcd$", base)
    if m:
        return int(m.group(1))
    return 10**9  # 无法解析的排后面

def main():
    if not os.path.exists(PDB_FILE):
        print("❌ 错误: 找不到 PDB 文件 {}".format(PDB_FILE))
        return
    if not os.path.isdir(TRAJ_DIR):
        print("❌ 错误: 找不到轨迹目录 {}".format(TRAJ_DIR))
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    # 1) 构建 residue pairs
    top = md.load(PDB_FILE).topology
    n_res = top.n_residues
    print("🧬 残基总数: {}".format(n_res))

    pairs = []
    for i in range(n_res):
        for j in range(i + 1 + GAP, n_res):
            pairs.append([i, j])
    pairs = np.array(pairs, dtype=np.int32)
    n_pairs = len(pairs)
    print("🔗 待计算残基对总数: {}".format(n_pairs))

    # 2) 搜索 DCD 文件（新命名）
    pat = os.path.join(TRAJ_DIR, "chain*_bulk_500ns.dcd")
    dcd_files = sorted(glob.glob(pat), key=natural_key_chain)

    if not dcd_files:
        print("❌ 没找到任何轨迹文件: {}".format(pat))
        return

    print("🔎 准备处理 {} 条轨迹：".format(len(dcd_files)))
    print("   e.g. {}".format(os.path.basename(dcd_files[0])))

    # 3) 统计 contact
    total_contact_counts = np.zeros(n_pairs, dtype=np.float64)
    total_frames = 0

    for dcd in tqdm(dcd_files, desc="Processing Trajectories", unit="traj"):
        # 分块读取轨迹
        for chunk in md.iterload(dcd, top=PDB_FILE, chunk=CHUNK, stride=STRIDE):
            # (frames, pairs)
            dists = md.compute_distances(chunk, pairs)
            contacts = (dists < RCUT)
            total_contact_counts += contacts.sum(axis=0)
            total_frames += chunk.n_frames

    if total_frames == 0:
        print("❌ 未读取到任何有效帧，请检查 stride/chunk 和轨迹文件。")
        return

    # 4) 归一化映射回矩阵
    avg_probs = total_contact_counts / float(total_frames)
    contact_map = np.zeros((n_res, n_res), dtype=np.float32)
    contact_map[pairs[:, 0], pairs[:, 1]] = avg_probs.astype(np.float32)
    contact_map[pairs[:, 1], pairs[:, 0]] = avg_probs.astype(np.float32)

    # 上三角 pair 的和（你定义的 sum_of_contacts）
    sum_of_contacts = float(np.sum(avg_probs))

    print("-" * 40)
    print("📊 统计分析完成 (基于 {} 帧, stride={}):".format(total_frames, STRIDE))
    print("   平均总接触数 (Sum of P_ij over i<j, gap>{}): {:.4f}".format(GAP, sum_of_contacts))
    print("-" * 40)

    # 5) 绘图（LogNorm）
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)

    plot_map = contact_map.copy().astype(np.float64)

    # mask 对角附近 gap
    for i in range(n_res):
        lower = max(0, i - GAP - 1)
        upper = min(n_res, i + GAP + 2)
        plot_map[i, lower:upper] = np.nan

    im = ax.imshow(
        plot_map,
        cmap="viridis",
        origin="lower",
        norm=LogNorm(vmin=VMIN, vmax=VMAX)
    )

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Contact Probability", fontsize=12)

    ax.set_title(
        "Global Ensemble Contact Map (260K)\nSum of Contacts = {:.1f}".format(sum_of_contacts),
        fontsize=15, pad=15
    )
    ax.set_xlabel("Residue Index", fontsize=12)
    ax.set_ylabel("Residue Index", fontsize=12)

    ticks = np.arange(0, n_res + 1, 40)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    out_img = os.path.join(OUT_DIR, "global_contact_log_bulk.png")
    out_npy = os.path.join(OUT_DIR, "global_contact_map_bulk.npy")

    plt.savefig(out_img, bbox_inches="tight")
    np.save(out_npy, contact_map)

    print("✅ 绘图已完成: {}".format(out_img))
    print("💾 数据已存为: {}".format(out_npy))

if __name__ == "__main__":
    main()
