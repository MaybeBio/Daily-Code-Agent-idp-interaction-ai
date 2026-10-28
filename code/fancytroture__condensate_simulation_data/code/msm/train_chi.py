#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import gc
import sys
import json
import numpy as np
import pandas as pd
# Select GPUs through CUDA_VISIBLE_DEVICES in the calling environment.
import tensorflow as tf

from tensorflow.keras.callbacks import ModelCheckpoint, CSVLogger

# 这两个线程设置在 tf.data 使用 from_tensor_slices 时影响不大，但保留无害
tf.config.threading.set_intra_op_parallelism_threads(8)
tf.config.threading.set_inter_op_parallelism_threads(8)

import model
from model import EarlyStopping 

# ==============================================================================
# WarmupAwareCheckpoint 类的定义 (保持不变)
# ==============================================================================
class WarmupAwareCheckpoint(ModelCheckpoint):
    def __init__(self, earlystop_callback: EarlyStopping, **kwargs):
        super().__init__(**kwargs)
        self.earlystop_callback = earlystop_callback

    def on_epoch_end(self, epoch, logs=None):
        current_epoch_num = epoch + 1
        start_epoch_num = self.earlystop_callback.start_epoch
        if current_epoch_num <= start_epoch_num:
            print(f" (Warmup Phase: Epoch {current_epoch_num}/{start_epoch_num}, Best Checkpoint Not Saved Yet)", end='')
            return
        if current_epoch_num == start_epoch_num + 1:
            self.best = self.earlystop_callback.best
            if np.isfinite(self.best):
                best_str = f"{self.best:.4f}"
            else:
                best_str = "Inf"
            print(f"\nWarmup finished. Checkpoint saving is now active. Initial best score: {best_str}")
        super().on_epoch_end(epoch, logs)

# ==============================================================================
# Eager Mode 检查 (保持不变)
# ==============================================================================
print(f"✅ TensorFlow v{tf.__version__}")
if not tf.executing_eagerly():
    try: tf.compat.v1.enable_eager_execution()
    except AttributeError: pass
print(f"✅ Eager Execution Enabled: {tf.executing_eagerly()}")

# ==============================================================================
# 参数与路径 (保持不变)
# ==============================================================================
attempt = int(sys.argv[1]) if len(sys.argv) > 1 else 0
n       = int(sys.argv[2]) if len(sys.argv) > 2 else 4
lag     = int(sys.argv[3]) if len(sys.argv) > 3 else 50
prefix = "combined"
width, depth = 256, 4
learning_rate, dropout, regularization = 1e-4, 0, 1e-5
dt = 0.01
train_ratio = 0.9
output_dir = "model_tdp43"
checkpoint_dir = "checkpoints_tdp43"
intermediate_dir = "intermediate_tdp43"
lag_dim = 10
os.makedirs(output_dir, exist_ok=True)
os.makedirs(checkpoint_dir, exist_ok=True)
os.makedirs(intermediate_dir, exist_ok=True)
model_path      = f"{output_dir}/model-ve-{prefix}-{n}-{attempt}-lag{lag}.h5"
checkpoint_best = f"{checkpoint_dir}/best-koop-{prefix}-n{n}-att{attempt}-lag{lag}.weights.h5"
checkpoint_last = f"{checkpoint_dir}/last-koop-{prefix}-n{n}-att{attempt}-lag{lag}.weights.h5"
log_path        = f"{checkpoint_dir}/log-{prefix}-n{n}-att{attempt}-lag{lag}.csv"
gen_save_path = f"{intermediate_dir}/model-idx-{prefix}-n{n}-{attempt}-lag{lag}.hdf5"
print(json.dumps({
    "attempt": attempt, "n": n, "prefix": prefix, "lag": lag, "dt": dt,
    "width": width, "depth": depth, "lr": learning_rate, "dropout": dropout,
    "regularization": regularization, "train_ratio": train_ratio
}, indent=2))

# ==============================================================================
# GPU 设置 (保持不变)
# ==============================================================================
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus: tf.config.experimental.set_memory_growth(gpu, True)
        print(f"✅ GPU memory growth = True, GPUs = {len(gpus)}")
    except Exception as e: print(f"⚠️ set_memory_growth failed:", e)
else: print("⚠️ 未检测到 GPU，将以 CPU 运行。")

# ==============================================================================
# 已完成检查 (简化)
# ==============================================================================
if os.path.exists(model_path):
    print(f"✅ 已存在最终模型，跳过训练。")
    sys.exit(0)

# ==============================================================================
# 5) 数据加载 (✅ *** 关键修改 ***)
# ==============================================================================
# train_chi.py

def load_trajs(feature_path, meta_path):
    if not os.path.exists(feature_path) or not os.path.exists(meta_path):
        sys.exit(f"❌ 缺少特征或元数据文件: {feature_path} 或 {meta_path}")
    
    # ------------------ 修正开始 ------------------
    print(f"🚀 正在将数据全部加载到 RAM 中 (这可能需要几秒钟)...")
    
    # ★★★ 你之前可能漏掉了下面这一行，或者变量名写错了 ★★★
    feats = np.load(feature_path) 
    
    print(f"✅ 特征文件加载完毕: shape={feats.shape}, dtype={feats.dtype}")
    # ------------------ 修正结束 ------------------

    # 下面是报错的第 111 行，现在 feats 已经定义了，就不会报错了
    sample_n = min(10000, feats.shape[0])
    sample = feats[:sample_n]
    print("🔎 DEBUG approx feats mean/std (first", sample_n, "frames):",
          float(sample.mean()), float(sample.std()))

    meta = np.load(meta_path)
    lengths = meta["lengths"]

    starts = np.concatenate([[0], np.cumsum(lengths[:-1])])
    trajs = [feats[start:start + L] for start, L in zip(starts, lengths)]
    return trajs


print("💾 加载全部轨迹...")
feature_path = f"{intermediate_dir}/{prefix}_tica_dim100_lag{lag_dim}_coords.npy"
meta_path    = f"{intermediate_dir}/{prefix}_meta.npz"
all_trajs = load_trajs(feature_path, meta_path)
np.random.seed(42)

# ==============================================================================
# DEBUG 模式 (保持不变)
# ==============================================================================
DEBUG_MODE = False
if DEBUG_MODE:
    print("\n" + "!"*56 + "\n!!!           警告：当前处于快速调试模式           !!!\n" + "!"*56 + "\n")
    import random
    trajs_for_training = random.sample(all_trajs, int(len(all_trajs) * 0.1) or 1)
    print(f"原始总轨迹数: {len(all_trajs)}, 采样后轨迹数: {len(trajs_for_training)}\n")
else:
    trajs_for_training = all_trajs

# ==============================================================================
# 6) 构建 KoopmanModel 和 DataGenerator
# ==============================================================================
print("🔧 构建 KoopmanModel 和 DataGenerator ...")
BATCH_SIZE = int(os.environ.get("BATCH_SIZE", "65536"))
print(f"✅ CHI batch_size={BATCH_SIZE}")
gen_all = model.DataGenerator(data=trajs_for_training, lag=lag, batch_size=BATCH_SIZE, ratio=train_ratio, dt=dt)
koop = model.KoopmanModel(n=n, network_lag=lag, verbose=1, nnargs=dict(width=width, depth=depth, learning_rate=learning_rate, regularization=regularization, dropout=dropout, batchnorm=True, lr_factor=5e-3))
koop.generator = gen_all
print("✅ KoopmanModel 和 DataGenerator 构建完成。")

# ==============================================================================
# 7) 断点续跑 (✅ *** 简化 ***)
# ==============================================================================
def safe_read_initial_epoch(csv_file):
    if not os.path.exists(csv_file): return 0, False
    try:
        df = pd.read_csv(csv_file)
        return len(df), True if not df.empty else False
    except Exception: return 0, False

initial_epoch, csv_append, loaded_weights = 0, False, False

if os.path.exists(gen_save_path):
    print(f"✅ 发现已保存的数据集索引，将加载以保证可复现性: {os.path.basename(gen_save_path)}")
    gen_all.load(gen_save_path)
else:
    print(f"✅ 未发现数据集索引，将创建新的划分并保存: {os.path.basename(gen_save_path)}")
    gen_all.save(gen_save_path)


for cp_path in [checkpoint_best, checkpoint_last]:
    if os.path.exists(cp_path):
        try:
            print(f"⏳ 尝试从 checkpoint 恢复权重: {os.path.basename(cp_path)}")
            koop.load_chi_weights(cp_path)
            loaded_weights = True
            print(f"✅ 成功恢复权重。")
            break
        except Exception as e: print(f"⚠️ 加载 {os.path.basename(cp_path)} 失败: {e}")

if loaded_weights:
    initial_epoch, csv_append = safe_read_initial_epoch(log_path)
    print(f"   - 将从 epoch {initial_epoch} 继续训练。")
elif os.path.exists(log_path):
    bak_path = f"{log_path}.bak.{pd.Timestamp.now().strftime('%Y%m%d%H%M%S')}"
    os.rename(log_path, bak_path)
    print(f"🧹 未加载到权重但发现旧日志，已重命名为: {os.path.basename(bak_path)}")

# ==============================================================================
# 8) 定义回调函数并训练 (保持不变)
# ==============================================================================
csv_logger = CSVLogger(log_path, append=csv_append)
early_stopping = EarlyStopping(monitor='val_metric_VAMP2', mode='max', min_delta=1e-3, patience=50, start_epoch=30, restore_best_weights=True, verbose=1)
checkpoint_best = WarmupAwareCheckpoint(earlystop_callback=early_stopping, filepath=checkpoint_best, monitor='val_metric_VAMP2', mode='max', save_best_only=True, save_weights_only=True, verbose=1)
checkpoint_last = ModelCheckpoint(filepath=checkpoint_last, save_best_only=False, save_weights_only=True, verbose=0)
callbacks = [checkpoint_best, checkpoint_last, csv_logger, early_stopping]

print(f"🚀 启动训练: n={n}, attempt={attempt}, initial_epoch={initial_epoch}")
koop.train_model(callbacks=callbacks, epochs=2000, initial_epoch=initial_epoch, verbose=1)
print("✅ chi 网络训练完成。")

# ==============================================================================
# 9) 保存模型 (保持不变)
# ==============================================================================
print("💾 保存最终模型 ...")
koop.save(model_path)
print(f"✅ 模型已保存: {model_path}")

del koop, all_trajs, gen_all
gc.collect()
print("🎉 脚本执行完毕。")
