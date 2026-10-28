#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, gc
# Select GPUs through CUDA_VISIBLE_DEVICES in the calling environment.
os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
import numpy as np
import model
from model import EarlyStopping

# ======================================================================
# 配置区（必须与 train_chi.py 严格一致）
# ======================================================================
attempt = int(sys.argv[1]) if len(sys.argv) > 1 else 0
n       = int(sys.argv[2]) if len(sys.argv) > 2 else 4
network_lag = int(sys.argv[3]) if len(sys.argv) > 3 else 200
prefix = "combined"
dt = 0.01
train_ratio = 0.9

# ⚠️ 必须与 train_chi.py 完全相同！
width = 256
depth = 4
learning_rate = 1e-4          # ← 关键：与 chi 一致！
dropout = 0
regularization = 1e-5
batchnorm = True

# ======================================================================
# 路径
# ======================================================================
intermediate_dir = "intermediate_tdp43"
ckpt_dir = "checkpoints_tdp43"
model_dir = "model_tdp43"
lag_dim = 10

feature_path = f"{intermediate_dir}/{prefix}_tica_dim100_lag{lag_dim}_coords.npy"
meta_path    = f"{intermediate_dir}/{prefix}_meta.npz"
gen_save_path = f"{intermediate_dir}/model-idx-{prefix}-n{n}-{attempt}-lag{network_lag}.hdf5"
chi_ckpt_path = f"{ckpt_dir}/best-koop-{prefix}-n{n}-att{attempt}-lag{network_lag}.weights.h5"
full_model_path = f"{model_dir}/model-ve-{prefix}-{n}-{attempt}-lag{network_lag}.h5"

for p in [feature_path, meta_path, gen_save_path, chi_ckpt_path]:
    if not os.path.exists(p):
        print(f"❌ 缺少文件: {p}")
        sys.exit(1)

print(f"\n🚀 开始训练辅助模型 (n={n}, lag={network_lag}, attempt={attempt})")

try:
    # 1. 创建模型
    koop = model.KoopmanModel(
        n=n, network_lag=network_lag, verbose=1,
        nnargs=dict(
            width=width, depth=depth, learning_rate=learning_rate,
            regularization=regularization, dropout=dropout,
            batchnorm=batchnorm, lr_factor=5e-3
        )
    )

    # 2. 加载数据
    all_frames = np.load(feature_path)
    meta = np.load(meta_path, allow_pickle=True)
    lengths = meta["lengths"]
    start_indices = np.concatenate([[0], np.cumsum(lengths[:-1])])
    trajectories = [all_frames[start:start+L] for start, L in zip(start_indices, lengths)]

    # 3. 初始化 DataGenerator（✅ 修复：补全参数）
    batch_size = int(os.environ.get("BATCH_SIZE", "65536"))
    print(f"✅ AUX batch_size={batch_size}")
    generator = model.DataGenerator(
        data=trajectories,
        lag=network_lag,          # ← 必须！
        batch_size=batch_size,    # ← 必须！
        ratio=train_ratio,
        dt=dt
    )
    generator.load(gen_save_path)
    koop.generator = generator

    # 4. 加载 chi 权重
    koop.load_chi_weights(chi_ckpt_path)
    print("✅ chi 权重加载成功")

    # 5. 定义早停回调（✅ 已添加）
    aux_early_stopping = EarlyStopping(
        monitor='val_loss',
        mode='min',
        min_delta=1e-3,
        patience=40,
        restore_best_weights=True,
        verbose=1
    )
    aux_callbacks = [aux_early_stopping]

    # 6. 训练辅助模型
    koop.train_auxiliary_model(callbacks=aux_callbacks)
    print("✅ 辅助模型训练完成")

    # 7. 保存完整模型
    koop.save(full_model_path)
    print(f"💾 模型已保存: {full_model_path}")

except Exception as e:
    print(f"\n❌ 错误: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
finally:
    gc.collect()

print("\n🎉 完成！")
