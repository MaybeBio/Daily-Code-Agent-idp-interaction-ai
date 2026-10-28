#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_npz_chain.py

从 extracted_500ns_bulk_chains/chainXXX_bulk_500ns.dcd 生成：
  - {out_prefix}_features.npy : (total_frames, n_pairs)  float32 memmap
  - {out_prefix}_meta.npz     : lengths, mean, std (float32)

假设：
  - 所有 chain*.dcd 使用同一个拓扑 PDB（默认 starts/seed_000.pdb）
  - PDB 中所有 CA 数一致且每个 DCD 的原子数与 PDB 匹配
"""

import os
import re
import glob
import argparse

import numpy as np
import mdtraj as md
from numpy.lib.format import open_memmap


def parse_args():
    ap = argparse.ArgumentParser(
        description="从 chainXXX_bulk_500ns.dcd 生成 combined_features.npy + combined_meta.npz"
    )
    ap.add_argument(
        "--pdb",
        default="starts/seed_000.pdb",
        help="拓扑 PDB（默认: starts/seed_000.pdb）",
    )
    ap.add_argument(
        "--traj-glob",
        default="extracted_500ns_bulk_chains/chain*_bulk_500ns.dcd",
        help="DCD glob（默认: extracted_500ns_bulk_chains/chain*_bulk_500ns.dcd）",
    )
    ap.add_argument(
        "--out-prefix",
        default="intermediate_tdp43/combined",
        help="输出前缀（默认: intermediate_tdp43/combined）",
    )
    ap.add_argument(
        "--skip-local",
        type=int,
        default=2,
        help="跳过 |i-j| <= skip_local 的残基对，默认 2",
    )
    ap.add_argument(
        "--chunk-frames",
        type=int,
        default=2000,
        help="分块计算距离时每块帧数（默认 2000）",
    )
    ap.add_argument(
        "--no-pbc",
        action="store_true",
        help="计算距离时不使用 PBC（periodic=False）",
    )
    return ap.parse_args()


def build_atom_pairs_ca(topology, skip_local=2):
    """
    从拓扑中找到所有 Cα 原子，并生成：
      - atom_pairs: 用于 md.compute_distances 的原子 index 对
      - res_pairs : 对应的残基 index 对（可选）
    条件：i_res < j_res 且 |i_res - j_res| > skip_local
    """
    ca_idx = topology.select("name CA")
    if ca_idx.size == 0:
        raise ValueError("拓扑中没有找到任何 Cα 原子 (select 'name CA')")

    res_idx = np.array([topology.atom(i).residue.index for i in ca_idx], dtype=int)

    atom_pairs = []
    res_pairs = []
    n_ca = len(ca_idx)
    for a in range(n_ca):
        for b in range(a + 1, n_ca):
            if abs(res_idx[b] - res_idx[a]) <= skip_local:
                continue
            atom_pairs.append((int(ca_idx[a]), int(ca_idx[b])))
            res_pairs.append((int(res_idx[a]), int(res_idx[b])))

    atom_pairs = np.array(atom_pairs, dtype=np.int32)
    res_pairs = np.array(res_pairs, dtype=np.int32)
    return atom_pairs, res_pairs


def main():
    args = parse_args()

    if not os.path.exists(args.pdb):
        raise RuntimeError("找不到 PDB：%s" % args.pdb)

    traj_files = sorted(glob.glob(args.traj_glob))
    if not traj_files:
        raise RuntimeError("traj_glob 没有匹配到任何 DCD：%r" % args.traj_glob)

    # 解析 chain id（仅用于打印/排序展示）
    pat_chain = re.compile(r"chain(\d+)")
    traj_info = []  # [(path, chain_id), ...]
    print("找到 DCD 轨迹（chain -> file）:")
    for f in traj_files:
        base = os.path.basename(f)
        m = pat_chain.search(base)
        if not m:
            print("  [警告] 文件名中未找到 'chainXXX'：%s，跳过" % f)
            continue
        chain_id = int(m.group(1))
        traj_info.append((f, chain_id))
        print("  chain %03d -> %s" % (chain_id, f))

    if not traj_info:
        raise RuntimeError("没有任何 DCD 通过命名解析出 chainXXX，请检查文件名模式。")

    # --- 参考拓扑构建 Cα 对 ---
    print("\n[参考拓扑] 使用 PDB:", args.pdb)
    t0 = md.load(args.pdb)
    top_ref = t0.topology
    ca_idx_ref = top_ref.select("name CA")

    atom_pairs, res_pairs = build_atom_pairs_ca(top_ref, skip_local=args.skip_local)
    n_pairs = atom_pairs.shape[0]
    print("Cα 原子数:", len(ca_idx_ref))
    print("特征对数 n_pairs:", n_pairs)
    print("前 5 个原子对 (atom index):", atom_pairs[:5])
    print("前 5 个残基对 (res index):", res_pairs[:5])

    # --- 第一轮：统计每条 DCD 的帧数 ---
    lengths = []
    total_frames = 0

    print("\n[第一轮] 统计每条 DCD 的帧数...")
    for f, chain_id in traj_info:
        traj = md.load(f, top=args.pdb)
        n_frames = traj.n_frames
        lengths.append(n_frames)
        total_frames += n_frames
        print(
            "  {:45s} chain={:03d}, n_frames={:7d}".format(
                os.path.basename(f), chain_id, n_frames
            )
        )
        del traj

    lengths = np.array(lengths, dtype=np.int64)
    print("总帧数 total_frames =", total_frames)

    # --- 创建特征 memmap 文件 ---
    out_features = args.out_prefix + "_features.npy"
    out_dir = os.path.dirname(out_features) or "."
    os.makedirs(out_dir, exist_ok=True)
    print("\n[第二步] 创建特征 memmap 文件:", out_features)
    X_mm = open_memmap(
        out_features,
        mode="w+",
        dtype="float32",
        shape=(total_frames, n_pairs),
    )

    # --- 第二轮：分块计算距离并写入 memmap ---
    use_pbc = not args.no_pbc
    sum_feat = np.zeros(n_pairs, dtype=np.float64)
    sum_sq_feat = np.zeros(n_pairs, dtype=np.float64)

    print("\n[第三步] 分块计算 Cα–Cα 距离并填充特征矩阵...")
    offset = 0
    CHUNK = int(args.chunk_frames)

    for (f, chain_id), n_frames in zip(traj_info, lengths):
        base = os.path.basename(f)
        print(
            "  处理 DCD: {:45s}  chain={:03d}, n_frames={:7d}".format(
                base, chain_id, n_frames
            )
        )
        traj = md.load(f, top=args.pdb)

        if traj.n_frames != n_frames:
            raise RuntimeError(
                "第二轮读取的帧数 ({}) 与第一轮记录的 n_frames ({}) 不一致，文件：{}".format(
                    traj.n_frames, n_frames, f
                )
            )

        for start in range(0, n_frames, CHUNK):
            end = min(n_frames, start + CHUNK)
            sub = traj[start:end]
            d_chunk = md.compute_distances(sub, atom_pairs, periodic=use_pbc)
            X_mm[offset + start : offset + end, :] = d_chunk.astype("float32")
            sum_feat += d_chunk.sum(axis=0)
            sum_sq_feat += (d_chunk ** 2).sum(axis=0)
            del sub, d_chunk

        offset += n_frames
        del traj

    if offset != total_frames:
        raise RuntimeError(
            "写入的总帧数 ({}) 与 total_frames ({}) 不一致！".format(offset, total_frames)
        )

    # --- 计算全局 mean / std ---
    print("\n[第四步] 计算 mean / std ...")
    mean = sum_feat / float(total_frames)
    var = sum_sq_feat / float(total_frames) - mean ** 2
    var = np.maximum(var, 1e-12)
    std = np.sqrt(var)

    # --- 保存 meta ---
    out_meta = args.out_prefix + "_meta.npz"
    print("\n[第五步] 保存 meta 到:", out_meta)
    np.savez(
        out_meta,
        lengths=lengths,
        mean=mean.astype("float32"),
        std=std.astype("float32"),
        # 如需保存配对信息可加：
        # res_pairs=res_pairs,
    )

    print("\n✅ 完成！")
    print("  特征矩阵:  {}  (shape = ({}, {}))".format(out_features, total_frames, n_pairs))
    print("  Meta 文件: {}".format(out_meta))


if __name__ == "__main__":
    main()
