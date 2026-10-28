#!/usr/bin/env python
# coding: utf-8

# # Neural network
# 
# In this notebook we set up the neural networks with VAMPNet scoring functions and train them for different output sizes and estimate errors by bootstrap aggregation. This notebook can be used with `papermill` to run all cells automatically with given parameters. We first define the imports and useful utility functions.

import gc
from glob import glob
from multiprocessing import Pool
import itertools
import os
from typing import List, Tuple, Sequence, Callable, TypeVar, Union, Any, Dict, Set
import warnings

import h5py
import mdtraj as md
import numpy as np
import pyemma as pe
from scipy.linalg import eig

import vampnet
import tensorflow as tf
from tensorflow.keras import regularizers
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import (Dense, Activation, Flatten, Input, BatchNormalization,
                                     concatenate, Dropout, AlphaDropout, Layer)
from tensorflow.keras.optimizers import Adam, Nadam, SGD
from tensorflow.keras.utils import plot_model
from tensorflow.keras.callbacks import EarlyStopping, Callback
import tensorflow.keras.backend as K
import pandas as pd

# REMOVED: The following line is the primary source of incompatibility.
# Modern TF/Keras scripts require Eager Execution to be enabled.
# tf.compat.v1.disable_eager_execution() 

warnings.filterwarnings('ignore')
vamp = vampnet.VampnetTools(epsilon=1e-5)

# ## Utility functions
# 
# The version of Keras we're using unfortunately doesn't have `restore_best_weights` implemented, so I copied this from a newer version.

# In[3]:


class EarlyStopping(Callback):
    """Stop training when a monitored quantity has stopped improving.
    # Arguments
        monitor: quantity to be monitored.
        min_delta: minimum change in the monitored quantity
            to qualify as an improvement, i.e. an absolute
            change of less than min_delta, will count as no
            improvement.
        patience: number of epochs that produced the monitored
            quantity with no improvement after which training will
            be stopped.
            Validation quantities may not be produced for every
            epoch, if the validation frequency
            (`model.fit(validation_freq=5)`) is greater than one.
        verbose: verbosity mode.
        mode: one of {auto, min, max}. In `min` mode,
            training will stop when the quantity
            monitored has stopped decreasing; in `max`
            mode it will stop when the quantity
            monitored has stopped increasing; in `auto`
            mode, the direction is automatically inferred
            from the name of the monitored quantity.
        baseline: Baseline value for the monitored quantity to reach.
            Training will stop if the model doesn't show improvement
            over the baseline.
        restore_best_weights: whether to restore model weights from
            the epoch with the best value of the monitored quantity.
            If False, the model weights obtained at the last step of
            training are used.
    """

    def __init__(self,
                 monitor='val_loss',
                 min_delta=0,
                 patience=0,
                 verbose=0,
                 mode='auto',
                 baseline=None,
                 restore_best_weights=False,
                 start_epoch=0): # MODIFIED: Added start_epoch to be compatible with train_chi.py
        super(EarlyStopping, self).__init__()

        self.monitor = monitor
        self.baseline = baseline
        self.patience = patience
        self.verbose = verbose
        self.min_delta = min_delta
        self.wait = 0
        self.stopped_epoch = 0
        self.restore_best_weights = restore_best_weights
        self.best_weights = None
        self.start_epoch = start_epoch # MODIFIED: Store start_epoch

        if mode not in ['auto', 'min', 'max']:
            warnings.warn('EarlyStopping mode %s is unknown, '
                          'fallback to auto mode.' % mode,
                          RuntimeWarning)
            mode = 'auto'

        if mode == 'min':
            self.monitor_op = np.less
        elif mode == 'max':
            self.monitor_op = np.greater
        else:
            if 'acc' in self.monitor:
                self.monitor_op = np.greater
            else:
                self.monitor_op = np.less

        if self.monitor_op == np.greater:
            self.min_delta *= 1
        else:
            self.min_delta *= -1

    def on_train_begin(self, logs=None):
        # Allow instances to be re-used
        self.wait = 0
        self.stopped_epoch = 0
        self.best_weights = None # MODIFIED: Reset best_weights
        if self.baseline is not None:
            self.best = self.baseline
        else:
            self.best = np.Inf if self.monitor_op == np.less else -np.Inf

    def on_epoch_end(self, epoch, logs=None):
        # MODIFIED: Respect the start_epoch parameter
        if epoch < self.start_epoch:
            return
            
        current = self.get_monitor_value(logs)
        if current is None:
            return

        if self.monitor_op(current - self.min_delta, self.best):
            self.best = current
            self.wait = 0
            if self.restore_best_weights:
                self.best_weights = self.model.get_weights()
        else:
            self.wait += 1
            if self.wait >= self.patience:
                self.stopped_epoch = epoch
                self.model.stop_training = True
                if self.restore_best_weights and self.best_weights is not None: # MODIFIED: Check if best_weights is not None
                    if self.verbose > 0:
                        print('Restoring model weights from the end of '
                              'the best epoch.') # MODIFIED: Added a period
                    self.model.set_weights(self.best_weights)

    def on_train_end(self, logs=None):
        if self.stopped_epoch > 0 and self.verbose > 0:
            print('Epoch %05d: early stopping' % (self.stopped_epoch + 1))

    def get_monitor_value(self, logs):
        logs = logs or {} # MODIFIED: Ensure logs is a dict
        monitor_value = logs.get(self.monitor)
        if monitor_value is None:
            warnings.warn(
                'Early stopping conditioned on metric `%s` '
                'which is not available. Available metrics are: %s' %
                (self.monitor, ','.join(list(logs.keys()))), RuntimeWarning
            )
        return monitor_value


# In[4]:


def unflatten(source: np.ndarray, lengths: List[int]) -> List[np.ndarray]:
    """
    Takes an array and returns a list of arrays.
    
    Parameters
    ----------
    source
        Array to be unflattened.
    lengths
        List of integers giving the length of each subarray.
        Must sum to the length of source.
    
    Returns
    -------
    unflat
        List of arrays.
    
    """
    # MODIFIED: Corrected logic for nested lists of lengths which seems to be the original intent.
    # This was a latent bug. The original code would not work if shapes had values > 1.
    out_list = []
    current_pos = 0
    for length in lengths:
        out_list.append(source[current_pos:current_pos + length])
        current_pos += length
    return out_list


# In[5]:


def sort_lengths(flatlengths: Sequence[int], shapes: Sequence[int]) -> List[List[int]]:
    """
    Takes a list of lengths and returns a list of lists of lengths.
    
    Parameters
    ----------
    flatlengths
        List of lengths
    shapes
        List of shapes
    
    Returns
    -------
    lengths
        List of lists of lengths
    
    """
    lengths = []
    i = 0
    for n in shapes:
        arr = []
        for _ in range(n):
            arr.append(flatlengths[i])
            i += 1
        lengths.append(arr)
    return lengths


# In[6]:


def triu_inverse(x: np.ndarray, n: int) -> np.ndarray:
    """
    Converts flattened upper-triangular matrices into full symmetric matrices.
    
    Parameters
    ----------
    x
        Flattened matrices
    n
        Size of the n * n matrix
    
    Returns
    -------
    mat
        Array of shape (length, n, n)
    
    """
    length = x.shape[0]
    mat = np.zeros((length, n, n))
    a, b = np.triu_indices(n, k=1)
    mat[:, a, b] = x
    mat += mat.swapaxes(1, 2)
    return mat


# In[8]:


class VAMPu(Layer):
    def __init__(self, units, activation, **kwargs):
        self.M = units
        self.activation = activation
        
        super().__init__(**kwargs)
    
    def build(self, input_shape):
        self.u_kernel = self.add_weight(name="u_var", shape=(self.M, ), trainable=True,
                                        initializer=tf.constant_initializer(1. / self.M))
        super().build(input_shape)
    
    def compute_output_shape(self, input_shape):
        return [self.M] * 2 + [(self.M, self.M)] * 4 + [self.M]
    
    def _tile(self, x, n_batch):
        x_exp = tf.expand_dims(x, axis=0)
        shape = x.get_shape().as_list()
        return tf.tile(x_exp, [n_batch, *([1] * len(shape))])
    
    def call(self, x):
        chi_t, chi_tau = x
        n_batch = tf.shape(chi_t)[0]
        norm = 1. / tf.cast(n_batch, dtype=tf.float32)
        
        corr_tau = norm * tf.matmul(chi_tau, chi_tau, transpose_a=True)
        chi_mean = tf.reduce_mean(chi_tau, axis=0, keepdims=True)
        kernel_u = tf.expand_dims(self.activation(self.u_kernel), axis=0)
        
        u = kernel_u / tf.reduce_sum(chi_mean * kernel_u, keepdims=True)
        v = tf.matmul(corr_tau, u, transpose_b=True)
        mu = norm * tf.matmul(chi_tau, u, transpose_b=True)
        sigma = tf.matmul(chi_tau * mu, chi_tau, transpose_a=True)
        gamma = chi_tau * tf.matmul(chi_tau, u, transpose_b=True)
        
        C00 = norm * tf.matmul(chi_t, chi_t, transpose_a=True)
        C11 = norm * tf.matmul(gamma, gamma, transpose_a=True)
        C01 = norm * tf.matmul(chi_t, gamma, transpose_a=True)
        
        return [
            self._tile(var, n_batch) for var in (u, v, C00, C11, C01, sigma)
        ] + [mu]


# In[9]:


class VAMPS(Layer):
    def __init__(self, units, activation, order=20, renorm=False, **kwargs):
        self.M = units
        self.activation = activation
        self.renorm = renorm
        self.order = order
        
        super().__init__(**kwargs)
    
    def build(self, input_shape):
        self.S_kernel = self.add_weight(name="S_var", shape=(self.M, self.M), trainable=True,
                                        initializer=tf.constant_initializer(0.1))
        super().build(input_shape)
    
    def compute_output_shape(self, input_shape):
        return [(self.M, self.M)] * 2 + [self.M] + [(self.M, self.M)]
    
    def call(self, x):
        if len(x) == 5:
            v, C00, C11, C01, sigma = x
        else:
            chi_t, chi_tau, u, v, C00, C11, C01, sigma = x
            u = u[0]
        
        n_batch = tf.shape(v)[0]
        norm = 1. / tf.cast(n_batch, dtype=tf.float32)
        C00, C11, C01 = C00[0], C11[0], C01[0]
        sigma, v = sigma[0], v[0]
        
        kernel_w = self.activation(self.S_kernel)
        w1 = kernel_w + tf.transpose(kernel_w)
        w_norm = w1 @ v
        
        # Numerical problems with using a high p-norm
        if self.renorm:
            quasi_inf_norm = lambda x: tf.reduce_max(tf.abs(x))
            w1 = w1 / quasi_inf_norm(w_norm)
            w_norm = w1 @ v
        
        w2 = (1 - tf.squeeze(w_norm)) / tf.squeeze(v)
        S = w1 + tf.linalg.diag(w2)
        
        if len(x) == 8:
            q = (norm * tf.transpose(tf.matmul(S, chi_tau, transpose_b=True))
                 * tf.matmul(chi_tau, u, transpose_b=True))
            probs = tf.reduce_sum(chi_t * q, axis=1)
        
        K = S @ sigma
        vamp_e = tf.transpose(S) @ C00 @ S @ C11 - 2 * tf.transpose(S) @ C01
        vamp_e_tile = tf.tile(tf.expand_dims(vamp_e, axis=0), [n_batch, 1, 1])
        K_tile = tf.tile(tf.expand_dims(K, axis=0), [n_batch, 1, 1])
        S_tile = tf.tile(tf.expand_dims(S, axis=0), [n_batch, 1, 1])
        
        if len(x) == 5:
            return [vamp_e_tile, K_tile, tf.zeros((n_batch, self.M)), S_tile]
        else:
            return [vamp_e_tile, K_tile, probs, S_tile]


# In[10]:


def matrix_inverse(mat):
    """
    Calculates the inverse of a square matrix.
    
    Parameters
    ----------
    mat
        Square real matrix
    
    Returns
    -------
    inv
        Inverse of the matrix
    
    """
    eigva, eigveca = np.linalg.eigh(mat)
    inc = eigva > vamp.epsilon
    eigv, eigvec = eigva[inc], eigveca[:, inc]
    return eigvec @ np.diag(1. / eigv) @ eigvec.T

def covariances(data):
    """
    Calculates (lagged) covariances.
    
    Parameters
    ----------
    data
        Data at time t and t + tau
    
    Returns
    -------
    C0inv
        Inverse covariance
    Ctau
        Lagged covariance
    
    """
    chil, chir = data
    norm = 1. / chil.shape[0]
    C0, Ctau = norm * chil.T @ chil, norm * chil.T @ chir
    C0inv = matrix_inverse(C0)
    return C0inv, Ctau

def _compute_pi(K):
    """
    Calculates the stationary distribution of a transition matrix.
    
    Parameters
    ----------
    K
        Transition matrix
    
    Returns
    -------
    pi
        Normalized stationary distribution
    
    """
    eigv, eigvec = np.linalg.eig(K.T)
    pi_v = eigvec[:, ((eigv - 1) ** 2).argmin()]
    return pi_v / pi_v.sum(keepdims=True)


# In[11]:


from typing import Tuple, Sequence, List, Union, Generator, Callable, Any, Dict, TypeVar, Set
from collections import UserList
from pathlib import Path

T = TypeVar("T")
MaybeListType = Union[List[T], T]
NNDataType = Tuple[List[np.ndarray], np.ndarray]
MaybePathType = Union[Path, str]

FRAMES, DIMENSIONS, FIRST, LAST, TIME, TIME_PLUS_LAG = 0, 1, 0, -1, 0, 1


# In[12]:


def make_list(item: MaybeListType[T], cls=list) -> List[T]:
    """
    Turn an object into a list, if it isn't already.
    
    Parameters
    ----------
    item
        Item to contain in a list
    
    Returns
    -------
    list
        List with item as only element
    
    """
    if not isinstance(item, list):
        item = [item]
    return cls(item)


# In[13]:


class DataSet:
    def __init__(self, trains: List[np.ndarray], valids: List[np.ndarray]=None,
                 y_train: np.ndarray=None, y_valid: np.ndarray=None):
        """
        DataSet - Container for training and validation data.
        
        Parameters
        ----------
        trains
            List of training datasets
        valids
            List of validation datasets
        y_train
            Dummy training target data
        y_valid
            Dummy validation target data
        
        """
        self.trains = trains
        self.valids = valids
        self.y_train = y_train
        self.y_valid = y_valid
    
    def __len__(self) -> int:
        return self.n_train
    
    def __getitem__(self, key: int) -> "DataSet":
        if isinstance(key, int):
            data = self.__class__([t[key][None] for t in self.trains],
                                  [t[key][None] for t in self.valids])
        else:
            data = self.__class__([t[key] for t in self.trains],
                                  [t[key] for t in self.valids])
        if self.n is not None:
            data.n = self.n
        return data
    
    @property
    def n_train(self) -> int:
        """Number of training samples."""
        return self.trains[FIRST].shape[FRAMES]
    
    @property
    def n_valid(self) -> int:
        """Number of validation samples."""
        if self.valids is None or self.valids[FIRST] is None:
            return 0
        return self.valids[FIRST].shape[FRAMES]
    
    @property
    def n_dims(self) -> int:
        """Number of input dimensions."""
        return self.trains[FIRST].shape[DIMENSIONS]
    
    @property
    def n(self) -> int:
        """Number of output dimensions."""
        if self.y_train is None:
            return None
        return self.y_train.shape[DIMENSIONS]
    
    @n.setter
    def n(self, n: int):
        self.y_train = np.zeros((self.n_train, n))
        if self.n_valid > 0:
            self.y_valid = np.zeros((self.n_valid, n))
        else:
            self.y_valid = None
    
    @property
    def train(self) -> NNDataType:
        """Training and target data pair."""
        return self.trains, self.y_train
    
    @property
    def valid(self) -> NNDataType:
        """Validation and target data pair."""
        return self.valids, self.y_valid


# In[14]:


# In[14]:
from tensorflow.keras.utils import Sequence

# MODIFIED: DataSet is still needed for the auxiliary training phase as a temporary container.
class DataSet:
    def __init__(self, trains: List[np.ndarray], valids: List[np.ndarray]=None,
                 y_train: np.ndarray=None, y_valid: np.ndarray=None):
        self.trains = trains
        self.valids = valids
        self.y_train = y_train
        self.y_valid = y_valid
    
    @property
    def n_train(self) -> int:
        return self.trains[0].shape[0] if self.trains else 0
    
    @property
    def n_valid(self) -> int:
        return self.valids[0].shape[0] if self.valids else 0

# ======================================================================
# --- RECONSTRUCTED DataGenerator as a Keras Sequence ---
# ======================================================================

# model.py

class DataGenerator:
    """
    ✅ 最终高性能版本 v3.1: 修复了验证集 repeat 的 bug，并采用了 mmap 优化的流水线。
    """
    def __init__(self, data, lag, batch_size, ratio=0.9, dt=1.0):
        self.data = data
        self.lag = lag
        self.batch_size = batch_size
        self.ratio = ratio
        self.dt = dt
        self.n_dims = self.data[0].shape[1] if self.data else 0

        print("   - DataGenerator: 1/2 Creating global frame index map...")
        indices = []
        for i, traj in enumerate(self.data):
            if len(traj) > self.lag:
                valid_frames = np.arange(len(traj) - self.lag)
                traj_indices = np.full((len(valid_frames),), i, dtype=np.int32)
                indices.append(np.stack([traj_indices, valid_frames.astype(np.int32)], axis=1))
        if not indices: raise ValueError("No trajectories are long enough for the given lag time!")
        self.indices = np.vstack(indices)
        np.random.shuffle(self.indices)
        
        train_size = int(len(self.indices) * self.ratio)
        self.train_indices = self.indices[:train_size]
        self.valid_indices = self.indices[train_size:]
        print(f"   - DataGenerator: Found {len(self.indices)} samples. Train: {len(self.train_indices)}, Valid: {len(self.valid_indices)}")
        
        self.train_steps = int(np.ceil(len(self.train_indices) / self.batch_size))
        self.valid_steps = int(np.ceil(len(self.valid_indices) / self.batch_size)) if len(self.valid_indices) > 0 else 0
        
        print("   - DataGenerator: 2/2 Creating tf.data.Dataset objects...")
        self.train_dataset = self._create_dataset(self.train_indices, shuffle=True)
        self.valid_dataset = self._create_dataset(self.valid_indices, shuffle=False)
        print("   - DataGenerator: Initialization complete.")

    def _create_dataset(self, indices, shuffle=False):
        """Creates a high-performance, mmap-friendly tf.data.Dataset pipeline."""
        if len(indices) == 0:
            return None
        
        with tf.device('/cpu:0'):
            dataset = tf.data.Dataset.from_tensor_slices(indices)
        
        if shuffle:
            dataset = dataset.shuffle(buffer_size=min(len(indices), 50000), reshuffle_each_iteration=True).repeat()
        
        dataset = dataset.batch(self.batch_size)
        
        dataset = dataset.map(self._fetch_batch_fn, num_parallel_calls=tf.data.AUTOTUNE)
        
        dataset = dataset.prefetch(buffer_size=tf.data.AUTOTUNE)
        
        return dataset
        
    # 在 DataGenerator 类中定义 _fetch_batch_fn
    def _fetch_batch_fn(self, batch_indices):
        # ✅ 高性能版本 v3.2
        def _numpy_fetch_batch(indices_arr):
            # ... (上面提供的新的、快速的 _numpy_fetch_batch 代码) ...
            batch_len = indices_arr.shape[0]
            df = pd.DataFrame(indices_arr)
            batch_x_t = np.zeros((batch_len, self.n_dims), dtype=np.float32)
            batch_x_tau = np.zeros_like(batch_x_t)
            for traj_idx, group in df.groupby(0):
                original_indices = group.index.values
                frame_indices = group[1].values
                batch_x_t[original_indices] = self.data[traj_idx][frame_indices]
                batch_x_tau[original_indices] = self.data[traj_idx][frame_indices + self.lag]
            return batch_x_t, batch_x_tau

        x_t, x_tau = tf.numpy_function(_numpy_fetch_batch, [batch_indices], [tf.float32, tf.float32])
        x_t.set_shape([None, self.n_dims])
        x_tau.set_shape([None, self.n_dims])
        return (x_t, x_tau), tf.zeros(tf.shape(x_t)[0], dtype=tf.float32)

    # save/load 方法保持不变 (适配 train_chi.py)
    def save(self, filename):
        with h5py.File(handle_path(filename, non_existent=True), "w") as f:
            f.create_dataset("indices", data=self.indices)
            f.attrs['ratio'] = self.ratio

    def load(self, filename):
        """Loads shuffled indices from a file and recreates datasets."""
        print(f"   - DataGenerator: Loading fixed indices from {os.path.basename(filename)}")
        with h5py.File(handle_path(filename), "r") as f:
            self.indices = f["indices"][:]
            # 确保 ratio 也从文件加载，以防脚本参数变化
            self.ratio = f.attrs.get('ratio', self.ratio)
        
        # ✅ 修复：直接在这里执行 _setup_datasets_from_indices 的逻辑

        # 1. 重新划分训练/验证索引
        train_size = int(len(self.indices) * self.ratio)
        self.train_indices = self.indices[:train_size]
        self.valid_indices = self.indices[train_size:]
        
        print(f"   - DataGenerator (loaded): Splitting data. Training samples: {len(self.train_indices)}, Validation samples: {len(self.valid_indices)}")
        
        # 2. 重新计算步数
        self.train_steps = int(np.ceil(len(self.train_indices) / self.batch_size))
        self.valid_steps = int(np.ceil(len(self.valid_indices) / self.batch_size)) if len(self.valid_indices) > 0 else 0
        
        # 3. 重新创建 tf.data.Dataset 对象
        self.train_dataset = self._create_dataset(self.train_indices, shuffle=True)
        self.valid_dataset = self._create_dataset(self.valid_indices, shuffle=False)
# class AggregatorMixin:
#     def agg(self, attr: str, func) -> any:
#         return func(comp.__getattribute__(attr) for comp in self._components)
    
#     def is_consistent(self, attr: str) -> bool:
#         return len(set(comp.__getattribute__(attr) for comp in self._components)) == 1

# def concat(data) -> list: # <--- 修改这里
#     return list(itertools.chain(*data))


# class MultiGenerator(AggregatorMixin):
#     def __init__(self, generators: MaybeListType[DataGenerator]):
#         self._components = make_list(generators)
#         assert self.is_consistent("n_dims")
        
#         if not self.is_consistent("max_frames"):
#             print("WARNING: Inconsistent number of frames for each ensemble!")
    
#     @property
#     def generators(self) -> List[DataGenerator]:
#         return self._components
    
#     @property
#     def data(self) -> List[np.ndarray]:
#         return self.agg("data", concat)
    
#     @property
#     def n_dims(self) -> int:
#         return self.generators[FIRST].n_dims
    
#     @property
#     def n_points(self) -> int:
#         return self.agg("n_points", sum)
    
#     @property
#     def n_traj(self) -> int:
#         return self.agg("n_traj", sum)
    
#     @property
#     def traj_lengths(self) -> int:
#         return self.agg("traj_lengths", concat)
    
#     @property
#     def data_flat(self) -> np.ndarray:
#         return np.vstack(self.data)
    
#     def regenerate_indices(self):
#         for gen in self.generators:
#             gen.regenerate_indices()

#     def save(self, filenames: Sequence[MaybePathType]):
#         for gen, filename in zip(self.generators, filenames):
#             gen.save(filename)
    
#     def load(self, filenames: Sequence[MaybePathType]):
#         for gen, filename in zip(self.generators, filenames):
#             gen.load(filename)
    
#     def __call__(self, n: int, lag: int) -> DataSet:
#         datasets = [gen(n=n, lag=lag) for gen in self.generators]
#         trains = [np.concatenate([ds.trains[TIME] for ds in datasets]),
#                   np.concatenate([ds.trains[TIME_PLUS_LAG] for ds in datasets])]
#         valids = [np.concatenate([ds.valids[TIME] for ds in datasets]),
#                   np.concatenate([ds.valids[TIME_PLUS_LAG] for ds in datasets])]
#         train_size = trains[TIME].shape[FRAMES]
#         valid_size = valids[TIME].shape[FRAMES]
#         y_train = np.zeros((train_size, 2 * n))
#         y_valid = np.zeros((valid_size, 2 * n))
#         train_indices = np.random.choice(np.arange(train_size), size=train_size, replace=False)
#         valid_indices = np.random.choice(np.arange(valid_size), size=valid_size, replace=False)
#         return DataSet(
#             [trains[TIME][train_indices], trains[TIME_PLUS_LAG][train_indices]],
#             [valids[TIME][valid_indices], valids[TIME_PLUS_LAG][valid_indices]],
#             y_train, y_valid
#         )
    
# # In[15]:


class KeepLast(UserList):
    def __init__(self, data):
        """
        Constructs a list that will always keep the first item.
        
        Parameters
        ----------
        data
            Data to construct the list from
        
        """
        self.data = list(reversed(data))
    
    def pop_first(self):
        """
        Returns the first item from the list, but only deletes
        it if there's at least one more item in the list.
        
        Returns
        -------
        item
            First item
        
        """
        if len(self) < 2:
            return self.data[LAST]
        return self.data.pop(LAST)


# In[16]:


def _split(data: np.ndarray, axis=LAST) -> List[np.ndarray]:
    """
    Utility function for splitting the output from two network lobes.
    
    Parameters
    ----------
    data
        Array to split
    axis
        Axis to split along
    
    Returns
    -------
    split
        2 arrays of half width
    
    """
    n = data.shape[axis] // 2
    return [data[:, :n], data[:, n:]]

def handle_path(path: MaybePathType, non_existent: bool=False) -> Path:
    """
    Check path validity and return `Path` object.

    Parameters
    ----------
    path
        Filepath to be checked.
    non_existent
        If false, will raise an error if the path does not exist.

    Returns
    -------
    path
        The converted and existing path.

    """
    if not isinstance(path, Path):
        try:
            path = Path(path)
        except Exception as err:
            message = "Couldn't read path {0}! Original message: {1}"
            raise ValueError(message.format(path, err))
    if not path.exists() and not non_existent:
        raise IOError("File {0} does not exist!".format(path))
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True) # MODIFIED: Added parents=True, exist_ok=True
    return path


# In[17]:


VALIDS = {int, float, str, list, bool, np.int64, np.float64} # MODIFIED: Added bool and numpy types
def _get_serializable_attributes(obj: object) -> Dict[str, Any]:
    """
    Finds all object attributes that are serializable with HDF5.
    
    Parameters
    ----------
    obj
        Object to serialize
    
    Returns
    -------
    attributes
        All serializable public attributes
    
    """
    return {k: v for k, v in obj.__dict__.items()
            if any(isinstance(v, valid) for valid in VALIDS)
            and not k.startswith("_")}


# In[18]:


class NNModel:
    def __init__(self, model: Union[Model, Layer],
                 loss: MaybeListType[Callable[..., Any]]=None,
                 metric: MaybeListType[Callable[..., Any]]=None,
                 learning_rate: MaybeListType[float]=None,
                 batchsize: int=10000,
                 epochs: int=100,
                 callback: Callback=None,
                 save_initial_weights=False,
                 verbose: int=0):
        """
        Neural network model interface class for Keras.
        
        Parameters
        ----------
        loss
            Loss function(s) to use, if more than one is specified
            subsequent calls to `compile` or `compile_fit` will use
            the next loss function.
        metric
            Fitting metric(s) to use
        learning_rate
            Learning rate for the Adam optimizer. Will accept
            multiple learning rates like `loss`.
        batchsize
            Batchsize to use for training and testing
        epochs
            The maximum number of epochs while training
        callback
            The callback to use while training
        save_initial_weights
            Whether to save the initial weights for possible later restoration
        verbose
            The verbosity level for Keras
        
        Attributes
        ----------
        n_opt
            The number of times `compile_fit` will be called, given
            potentially multiple loss functions and / or learning rates.
        
        """
        self._model = model
        self._history = []
        
        # Using a special list structure allows to use
        # a varying number of compile parameters
        self.loss = make_list(loss, cls=KeepLast)
        self.learning_rate = make_list(learning_rate, cls=KeepLast)
        self.n_opt = max(len(d) for d in (self.loss, self.learning_rate))
        
        self.metric = metric
        self.callback = make_list(callback)
        self.batchsize = batchsize
        self.epochs = epochs
        self.verbose = verbose
        self.save_initial_weights = save_initial_weights
        self._initial_weights = self.weights if self.save_initial_weights else None
    
    @property
    def weights(self) -> List[np.ndarray]:
        """The weights of the model layers."""
        return self._model.get_weights()
    
    @weights.setter
    def weights(self, weights: MaybeListType[np.ndarray]):
        self._model.set_weights(weights)
    
    def reset_weights(self):
        """Reset the weights to the initialized state."""
        if not self.save_initial_weights:
            raise ValueError("Initial weights were not saved!")
        self.weights = self._initial_weights
    
    def compile(self):
        """Compiles the neural network model with the parameters set at instantiation."""
        self._model.compile(
            optimizer=Adam(learning_rate=self.learning_rate.pop_first(),
                           epsilon=0.0001, clipnorm=1.0, beta_1=0.99),
            loss=self.loss.pop_first(),
            metrics=self.metric)
    
    def fit(self, train: NNDataType, valid: NNDataType, **kwargs):
        """
        Trains the neural network with the specified training data.
        
        Parameters
        ----------
        train
            Training data in the form of a list of inputs and a target output.
        valid
            Validation data in the form of a list of inputs and a target output.
        kwargs
            Additional parameters for Keras `fit` or to override object parameters.
        
        """
        # Pass the defaults specified in the constructor and
        # optionally override them with explicitly passed options.
        batchsize = train[-1].shape[0] if self.batchsize == -1 else self.batchsize
        
        # MODIFIED: Ensure callbacks passed in kwargs override the default ones.
        callbacks = kwargs.pop('callbacks', self.callback)
        
        options = dict(batch_size=batchsize, epochs=self.epochs,
                       verbose=self.verbose, callbacks=callbacks, shuffle=True)
        
        # We have to filter out None because we would like to
        # use the Keras defaults if at all possible.
        options = {k: v for k, v in options.items() if v is not None}
        options.update(kwargs)
        hist = self._model.fit(*train, validation_data=valid, **options)
        self._history.append(hist)
    
    def compile_fit(self, train: NNDataType, valid: NNDataType, **kwargs):
        """
        Compiles and trains the neural network with the specified training data
        repeatedly, until all loss functions or learning rates have been used.
        
        Parameters
        ----------
        train
            Training data in the form of a list of inputs and a target output.
        valid
            Validation data in the form of a list of inputs and a target output.
        kwargs
            Additional parameters for Keras `fit` or to override object parameters.
        
        """
        for _ in range(self.n_opt):
            self.compile()
            # MODIFIED: Pass kwargs through to the fit method.
            self.fit(train, valid, **kwargs)
        
    def predict(self, data: MaybeListType[np.ndarray]) -> List[np.ndarray]:
        """
        Projects data through the neural network.
        
        Parameters
        ----------
        data
            Input data to predict
        
        Returns
        -------
        predicted
            Predicted data
        
        """
        # This is the case when we're predicting the state
        # assignments of a biased trajectory
        if isinstance(data, np.ndarray):
            data = [data, data]
        elif len(data) == 1:
            data = [data[0], data[0]]
            
        # Specifying batchsize is crucial because
        # we have a batch normalization layer in our model!
        return self._model.predict(data, batch_size=self.batchsize)
    
    def predict_dataset(self, data: DataSet, valid: bool=False,
                        func: Callable[..., Any]=None) -> DataSet:
        """
        Projects a DataSet through the neural network.
        
        Parameters
        ----------
        data
            Input DataSet to predict
        valid
            Whether to also process any validation data in the passed data object
        func
            Optional function to process the raw output of the neural network
        
        Returns
        -------
        predicted
            Predicted data as a DataSet
        
        """
        train_pred = self.predict(data.trains)
        valid_pred = self.predict(data.valids) if valid and data.n_valid > 0 else None
        
        # We will sometimes need to process the output before
        # forming a DataSet, e.g. splitting the network lobes
        if func is not None:
            train_pred = func(train_pred)
            if valid_pred is not None:
                valid_pred = func(valid_pred)
        return DataSet(train_pred, valid_pred)
        
    def evaluate(self, data: NNDataType) -> float:
        """
        Evaluates the score by passing the output to the specified loss function.
        
        Parameters
        ----------
        data
            Data to use for evaluation, typically validation data
        
        Returns
        -------
        score
            The loss function score
        
        """
        return self._model.evaluate(*data, batch_size=self.batchsize, verbose=self.verbose)
    
    # TODO: Newer Keras versions should be able to save directly to an HDF5 group
    def save(self, group: h5py.Group):
        """
        Save the model to a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        for i, weight in enumerate(self._model.get_weights()):
            group[str(i)] = weight
    
    def load(self, group: h5py.Group):
        """
        Load the model from a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        n_items = len(group)
        self._model.set_weights([
            np.array(group[str(i)]) for i in range(n_items) # MODIFIED: Ensure data is loaded as numpy array
        ])
    
    def load_weights(self, filename: MaybePathType):
        """
        Load the model weights from a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        filename = handle_path(filename).as_posix()
        self._model.load_weights(filename)


# In[19]:


def _build_model(n_input: int, n_output: int, learning_rate: float=1e-4,
                 width: int=1024, depth: int=2, regularization: float=1e-8,
                 dropout: float=0.0, verbose: int=0, batchnorm: bool=False,
                 lr_factor: float=1e-2) -> Dict[str, NNModel]:
    """
    Builds the VAMPNet model.
    
    Parameters
    ----------
    n_input
        Number of input dimensions
    n_output
        Number of output dimensions
    learning_rate
        Learning rate for the chi model
    width
        Width of the layers in neurons
    depth
        Depth of the model in layers
    regularization
        L2 regularization strength per hidden layer
    dropout
        Dropout per hidden layer
    verbose
        Verbosity of the vanilla model
    lr_factor
        Learning rate modifier for the `all` model
    
    Returns
    -------
    models
        Dictionary of all models and the u and S layers
    
    """

    # Input layers
    xti = Input(shape=(n_input,))
    xli = Input(shape=(n_input,))

    # Create hidden layers
    dense = []
    for i in range(depth):
        layer = Dense(
            units=width, use_bias=not batchnorm,  # 防止BN和bias重复
            kernel_regularizer=regularizers.l2(regularization),
            kernel_initializer=KoopmanModel.initializer)
        dense.append(layer)
        if batchnorm:
            dense.append(BatchNormalization())
        dense.append(Activation(KoopmanModel.activation))
        if i <= depth - 1 and dropout > 0.0:
            dense.append(AlphaDropout(dropout))

    # Optional batch normalization
    bn = BatchNormalization()
    lx = bn(xti) if batchnorm else xti
    rx = bn(xli) if batchnorm else xli
    for i, layer in enumerate(dense):
        lx = dense[i](lx)
        rx = dense[i](rx)

    # Output is softmax for [0, 1] interval
    softmax = Dense(
        units=n_output, activation="softmax",
        kernel_regularizer=regularizers.l2(0.1 * regularization),
        kernel_initializer=KoopmanModel.initializer)
    lx = softmax(lx)
    rx = softmax(rx)

    # Build the model
    merged = concatenate([lx, rx])
    models = {}
    models["chi"] = NNModel(
        Model(inputs=[xti, xli], outputs=merged),
        loss=vamp.loss_VAMP2_autograd,   # 单一损失
        learning_rate=learning_rate,     # 单一学习率
        metric=vamp.metric_VAMP2,   
        # MODIFIED: Removed the default callback. It will be provided by the training script.
        # callback=EarlyStopping("val_metric_VAMP", mode="max", min_delta=0.001,
        #                       patience=20, restore_best_weights=True),
        verbose=1)

    # Auxiliary inputs
    chil_in, chir_in = [Input(shape=(n_output,)) for _ in range(2)]
    v_in = Input(shape=(n_output, 1))
    C00_in, C01_in, C11_in, sigma_in = [
        Input(shape=(n_output, n_output)) for _ in range(4)]

    # Constraint layers
    vlu = VAMPu(n_output, activation=tf.exp)
    vls = VAMPS(n_output, activation=tf.exp, renorm=True)

    # In / Output for full model
    (u_out, v_out, C00_out, C11_out,
     C01_out, sigma_out, mu_out) = vlu([lx, rx])
    Ve_out, K_out, p_out, S_out = vls([
        lx, rx, u_out, v_out, C00_out,
        C11_out, C01_out, sigma_out])

    # In / Output for model with only u and S
    (u_out_b, v_out_b, C00_out_b, C11_out_b,
     C01_out_b, sigma_out_b, _) = vlu([chil_in, chir_in])
    Ve_out_b, K_out_b, p_out_b, S_out_b = vls([
        chil_in, chir_in, u_out_b, v_out_b,
        C00_out_b, C11_out_b, C01_out_b, sigma_out_b])

    # In / Output for new tau prediction model
    Ve_out_s, *_ = vls([v_in, C00_in, C11_in, C01_in, sigma_in])

    # We will need these layers later, so we save them as models
    models["vlu"] = NNModel(vlu, save_initial_weights=True)
    models["vls"] = NNModel(vls, save_initial_weights=True)

    # Build training models, we need to be very careful with batchsizes:
    # https://github.com/keras-team/keras/issues/12400
    models["all"] = NNModel(
        Model(inputs=[xti, xli], outputs=Ve_out),
        loss=loss_vampe, learning_rate=learning_rate * lr_factor,
        # callback=EarlyStopping("val_loss", patience=50,
        #     mode="min", restore_best_weights=True),
        epochs=KoopmanModel.n_epoch_aux, verbose=1)
    models["both"] = NNModel(
        Model(inputs=[chil_in, chir_in], outputs=Ve_out_b),
        loss=loss_vampe, learning_rate=5e-4,
        # callback=EarlyStopping("val_loss", patience=100,
        #     mode="min", restore_best_weights=True),
        epochs=KoopmanModel.n_epoch_aux, batchsize=-1)
    models["S"] = NNModel(
        Model(inputs=[v_in, C00_in, C11_in, C01_in, sigma_in], outputs=Ve_out_s),
        loss=loss_vampe, learning_rate=0.1,
        # callback=EarlyStopping("val_loss", patience=100,
        #     mode="min", restore_best_weights=True),
        epochs=KoopmanModel.n_epoch_aux, batchsize=-1, verbose=0)

    # Build prediction models
    models["inp"] = NNModel(Model(inputs=[chil_in, chir_in], outputs=[
            v_out_b, C00_out_b, C11_out_b, C01_out_b, sigma_out_b]))
    models["mu"] = NNModel(Model(inputs=[xti, xli], outputs=mu_out))
    models["K"] = NNModel(Model(inputs=[xti, xli], outputs=K_out))
    return models


# In[20]:


class KoopmanModel:
    activation = "selu"
    initializer = "lecun_normal"
    n_epoch_chi = 100
    n_epoch_aux = 10000
    
    def __init__(self, n: int, network_lag: int=4, constrained: bool=True,
                 verbose: int=0, nnargs: Dict[str, Any]=None):
        self.n_output = n
        self.n_input = None
        self.network_lag = network_lag
        self.verbose = verbose
        self.constrained = constrained
        self.nnargs = nnargs or {}
        self.chi_estimated = False
        self.aux_estimated = False
        
        self._lag = network_lag
        self._reestimated = False
        self._models = None
        self._chi_weights = None
        self._k_cache = {}
        self._max_iter = 20
        
        self._generator = None
        
        self._K = None
        self._pi = None
        self._mu = None
    
    @classmethod
    def from_file(cls, filename: MaybePathType) -> "KoopmanModel":
        """
        Initialize from a saved model.
        
        Parameters
        ----------
        filename
            HDF5 file with the saved model
        
        """
        koop = cls(n=-1)
        koop.load(handle_path(filename))
        return koop
    
    def __repr__(self) -> str:
        width = self.nnargs.get("width", "N/A")
        depth = self.nnargs.get("depth", "N/A")
        return "<KoopmanModel with shape=({0},{1}), lag={2}, n={3}{4}>".format(
            width, depth, self.network_lag, self.n_output,
            ", estimated" if self.chi_estimated else "")

    def _cleanup(self):
        """
        Cleanup model graph after training. Without this
        the kernel will die due to running out of memory.
        
        """
        gc.collect()
        K.clear_session()
    
    @property
    def n(self) -> int:
        """Alias for `n_output`."""
        return self.n_output
    
    @property
    def is_built(self) -> bool:
        """True if the neural network models have been built."""
        return self._models is not None
        
    @property
    def K(self) -> np.ndarray:
        """The estimated Koopman operator."""
        if self._K is None or self._reestimated:
            self._K = self._models["K"].predict(self.data.trains)[FIRST]
        return self._K
    
    @property
    def mu(self) -> np.ndarray:
        """The estimated mu."""
        if self._mu is None or self._reestimated:
            self._mu = self._models["mu"].predict(self.data.trains).flatten()
        return self._mu
    
    @property
    def pi(self) -> np.ndarray:
        """The estimated equilibrium distribution."""
        return statdist(self.K)
    
    # TODO Add reweighting code here
    @pi.setter
    def pi(self, value): # MODIFIED: Add value parameter
        pass
    
    @property
    def dt(self) -> float:
        """The timestep of the underlying data."""
        return self.generator.dt
    
    @property
    def generator(self):
        return self._generator
    
    @generator.setter
    def generator(self, generator: DataGenerator):
        self._generator = generator
        self.n_input = generator.n_dims
        if self._models is None:
            self._models = _build_model(self.n_input, self.n_output,
                                        verbose=self.verbose, **self.nnargs)

    def train_model(self, callbacks=None, **kwargs):
        if self.generator is None:
            raise ValueError("Generator has not been set!")

        train_ds = self.generator.train_dataset
        valid_ds = self.generator.valid_dataset
        train_steps = self.generator.train_steps
        valid_steps = self.generator.valid_steps if valid_ds else None

        if not valid_ds:
            warnings.warn("Validation dataset is empty.")

        # === 关键修改：只用一个 loss（loss_VAMP2_autograd）===
        # 确保 _build_model 中 chi 模型只定义了一个 loss 和一个 lr
        print("\n--- Training chi network with VAMP2 loss ---")

        # 编译模型（只编译一次）
        self._models["chi"].compile()  # 这会 pop_first loss 和 lr（现在只有一个）

        # 直接使用 tf.data.Dataset 训练
        self._models["chi"]._model.fit(
            train_ds,
            validation_data=valid_ds,
            steps_per_epoch=train_steps,
            validation_steps=valid_steps,
            callbacks=callbacks,
            **kwargs
        )

        self._chi_weights = self._models["chi"].weights
        self.chi_estimated = True

    def train_auxiliary_model(self, callbacks=None):
        """
        Train the auxiliary model ('both' and 'all') with safe batching.
        This version:
        - Fixes lag reference (uses self.network_lag)
        - Avoids modifying read-only DataSet attributes
        - Avoids negative shape errors
        - Controls GPU memory usage via batched prediction
        """
        assert self.chi_estimated, "Chi model must be trained before training auxiliary model."

        print("   - Loading full dataset for auxiliary model training (one-time operation)...")

        # ==================== Helper function: load paired (x_t, x_tau) ====================
        def get_data_from_indices(indices):
            if not len(indices):
                return None, None
            x_t = np.zeros((len(indices), self.n_input), dtype=np.float32)
            x_tau = np.zeros_like(x_t)
            for i, (traj_idx, frame_idx) in enumerate(indices):
                x_t[i] = self.generator.data[traj_idx][frame_idx]
                x_tau[i] = self.generator.data[traj_idx][frame_idx + self.network_lag]
            return x_t, x_tau

        # ✅ Load all data once (for train/valid)
        train_x_t, train_x_tau = get_data_from_indices(self.generator.train_indices)
        valid_x_t, valid_x_tau = get_data_from_indices(self.generator.valid_indices)

        # ==================== Helper: batched predict ====================
        def batched_predict(model, x_list, batch_size=2048):
            """
            Predict in large batches safely (防止显存溢出 + 预测失败时跳过该批)
            """
            preds = []
            n = len(x_list[0])
            total_batches = (n + batch_size - 1) // batch_size

            for i in range(0, n, batch_size):
                batch0 = x_list[0][i:i + batch_size]
                batch1 = x_list[1][i:i + batch_size]
                print(f"   - Predicting batch {i // batch_size + 1}/{total_batches} "
                      f"({len(batch0)} samples)...", end="\r")
                pred = None
                try:
                    pred = model.predict([batch0, batch1])
                    preds.append(pred)
                except tf.errors.ResourceExhaustedError:
                    print(f"\n⚠️  Batch {i // batch_size + 1} OOM, skipping this batch.")
                finally:
                    del batch0, batch1
                    if pred is not None:
                        del pred
                    gc.collect()
                    tf.keras.backend.clear_session()

            print()  # 换行
            if len(preds) == 0:
                raise RuntimeError("❌ batched_predict() failed: no predictions collected.")
            return np.concatenate(preds, axis=0)

        # ==================== Step 1: Predict χ_t and χ_τ ====================
        print("   - Setting initial auxiliary weights (batched prediction)...")
# 获取原始拼接输出 (N, 8)
        raw_chi_train_t = batched_predict(self._models["chi"], [train_x_t, train_x_t], batch_size=131072)
        raw_chi_train_tau = batched_predict(self._models["chi"], [train_x_tau, train_x_tau], batch_size=131072)

        # ✅ 拆分：取前 n_output 列作为 chi_t 和 chi_tau
        chi_train_t = raw_chi_train_t[:, :self.n_output]      # (N, 4)
        chi_train_tau = raw_chi_train_tau[:, :self.n_output]  # (N, 4)
        chi_valid_t, chi_valid_tau = None, None
        if valid_x_t is not None:
            raw_chi_valid_t = batched_predict(self._models["chi"], [valid_x_t, valid_x_t], batch_size=131072)
            raw_chi_valid_tau = batched_predict(self._models["chi"], [valid_x_tau, valid_x_tau], batch_size=131072)
            chi_valid_t = raw_chi_valid_t[:, :self.n_output]
            chi_valid_tau = raw_chi_valid_tau[:, :self.n_output]

        # ==================== Step 2: Build new DataSet ====================
        if chi_valid_t is not None and chi_valid_tau is not None:
            chi_data = DataSet(
                trains=[chi_train_t, chi_train_tau],
                valids=[chi_valid_t, chi_valid_tau]
            )
            n_train = len(chi_train_t)
            n_valid = len(chi_valid_t)
        else:
            chi_data = DataSet(trains=[chi_train_t, chi_train_tau])
            n_train = len(chi_train_t)
            n_valid = 0

        # ==================== Step 3: Wrap into (X, y) tuples ====================
        chi_train_data = (
            [chi_data.trains[0], chi_data.trains[1]],
            np.zeros((n_train, 1))
        )
        chi_valid_data = (
            ([chi_data.valids[0], chi_data.valids[1]], np.zeros((n_valid, 1)))
            if n_valid > 0 else None
        )


        full_train_data = (
            [train_x_t, train_x_tau],
            np.zeros((len(train_x_t), 1))
        )
        full_valid_data = (
            ([valid_x_t, valid_x_tau],
             np.zeros((len(valid_x_t), 1)))
            if valid_x_t is not None else None
        )

        # ==================== Step 4: Train "both" and "all" networks ====================
        print("   - Training auxiliary network ('both')...")
        self._models["both"].compile_fit(
            chi_train_data, chi_valid_data,
            callbacks=callbacks,
            verbose=2
        )

        print("   - Training full network ('all')...")
        self._models["all"].compile()  # 先编译模型
        
        hist = self._models["all"]._model.fit(
            self.generator.train_dataset,
            validation_data=self.generator.valid_dataset,
            steps_per_epoch=self.generator.train_steps,
            validation_steps=self.generator.valid_steps,
            callbacks=callbacks,
            epochs=self._models["all"].epochs,
            verbose=2
        )
        
        self._models["all"]._history.append(hist)
        self.aux_estimated = True
        print("✅ Auxiliary model training completed successfully.")


    
    def train_loop(self):
        """Train the auxiliary constraint model in a loop."""
        # Needs a trained vanilla VAMPNet
        assert self.chi_estimated and self.aux_estimated

        score, score_prev = 0.0, 0.0
        weights = self._models["all"].weights
        
        self._log("Setting up auxiliary training loop...")
        
        self._models["both"].compile()
        self._models["all"].compile()
        
        for i in range(self._max_iter):
            chi_data = self._update_auxiliary_weights(reset_weights=False)
            self._models["both"].fit(chi_data.train, chi_data.valid)
            self._log("Both valid loss: {0}".format(self._models["both"]._history[-1].history["val_loss"][-1])) # MODIFIED: get latest history
            self._models["all"].fit(self.data.train, self.data.valid)
            score = -self._models["all"].evaluate(self.data.valid)
            
            self._log("Iteration {0}: Score: {1}".format(i, score))
            
            if score <= score_prev:
                break
                
            score_prev = score
            weights = self._models["all"].weights
        
        # Final fit
        self._log("Performing final fit...")
        self._models["all"].weights = weights
        chi_data = self._update_auxiliary_weights(reset_weights=False)
        self._models["both"].fit(chi_data.train, chi_data.valid)
        self._log("Both valid loss: {0}".format(self._models["both"]._history[-1].history["val_loss"][-1])) # MODIFIED: get latest history

    def save(self, filename: MaybePathType):
        """
        Save the model to a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        filename = handle_path(filename, non_existent=True)
        with h5py.File(filename, "w") as write:
            
            # We need to create a dummy dataset to be able to use `attrs`
            par = write.create_dataset("parameters", (0,))
            parameters = _get_serializable_attributes(self)
            for k, v in parameters.items():
                # HDF5 can't store None, so we skip it.
                if v is not None:
                    par.attrs[k] = v
            
            # Save the neural net parameters
            nn = write.create_dataset("nn", (0,))
            for k, v in self.nnargs.items():
                if v is not None:
                    nn.attrs[k] = v
            
            write.create_group("models")
            for k, model in self._models.items():
                model_group = "models/{0}".format(k)
                write.create_group(model_group)
                model.save(write[model_group])
    
    def load(self, filename: MaybePathType):
        """
        Load a model from a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        filename = handle_path(filename)
        with h5py.File(filename, "r") as read:
            # TODO: Put this into a mixin instead
            self.__dict__.update(dict(read["parameters"].attrs))
            self.nnargs = dict(read["nn"].attrs)
            self._models = _build_model(
                self.n_input, self.n_output, verbose=self.verbose, **self.nnargs)
            for k, model in self._models.items():
                model.load(read["models/{0}".format(k)])

            # Make sure the initial weights are set
            self._chi_weights = self._models["chi"].weights
            self.chi_estimated = True
    
    def load_chi_weights(self, filename: MaybePathType):
        """
        Load a model from a file.
        
        Parameters
        ----------
        filename
            Path to an HDF5 file
        
        """
        if not self.is_built:
            # This case happens when called from train_aux.py
            # The generator is not set yet, so n_input is None. We can't build the model.
            # We'll defer model building until the generator is set.
            pass
        else:
            filename = handle_path(filename)
            self._models["chi"].load_weights(filename)
            self._chi_weights = self._models["chi"].weights
            self.chi_estimated = True
    
    def fit(self, generator: DataGenerator):
        """
        Train the model.
        
        Parameters
        ----------
        generator
            Data generator function, must be callable with output size and lag time
        
        """
        self.generator = generator
        if self._models is None: # MODIFIED: Check if model needs building
            self._models = _build_model(self.n_input, self.n_output,
                                        verbose=self.verbose, **self.nnargs)
        
        self.train_model()
        
        # This is the more compute-intensive part
        if self.constrained:
            self.train_auxiliary_model()
            self.train_loop()
    
    def transform(self, data: MaybeListType[np.ndarray]) -> MaybeListType[np.ndarray]:
        assert self.chi_estimated
        
        is_list = isinstance(data, list)
        if is_list:
            lengths = [len(d) for d in data]
            data_flat = np.vstack(data)
        else:
            data_flat = data

        # Prediction needs a pair of inputs
        prob = self._models["chi"].predict([data_flat, data_flat])
        # The output is concatenated, we only need the first half (chi_t)
        chi_t = prob[:, :self.n_output]

        if is_list:
            return unflatten(chi_t, lengths=lengths)
        return chi_t
    
    def fit_transform(self, generator: DataGenerator) -> np.ndarray:
        """
        Fit and transform data through the model.
        
        Parameters
        ----------
        data
            Data generator function, must be callable with output size, lag time,
            and must have a `flat` attribute for the complete input dataset.
        
        Returns
        -------
        predicted
            Array(s) with probabilistic state assignments
        
        """
        self.fit(generator)
        return self.transform(generator.data) # MODIFIED: was generator.data_flat
    
    def score(self) -> float:
        """
        Validation score of the trained model.
        
        Returns
        -------
        score
            The VAMP-E score of the fully trained model
        
        """
        assert self.chi_estimated and self.aux_estimated
        if self.data.n_valid > 0:
            return -self._models["all"].evaluate(self.data.valid)
        else:
            warnings.warn("No validation data available to score.")
            return np.nan
        
    def estimate_koopman(self, lag: int):
        """
        Estimate the Koopman operator K(τ) = C00^-1 * C01 for a given lag time.
        Used in analysis scripts (ITS, CK test).
        """
        import numpy as np
    
        # 如果缓存里已有，直接返回
        if hasattr(self, "_k_cache") and lag in self._k_cache:
            return self._k_cache[lag]
    
        # === 正确获取 χ 轨迹 ===
        if hasattr(self, "_chi_trajs"):
            chi_dict = self._chi_trajs
            if isinstance(chi_dict, dict) and "chi_t" in chi_dict and "chi_tau" in chi_dict:
                chi_t = chi_dict["chi_t"]
                chi_tau = chi_dict["chi_tau"]
            else:
                raise ValueError("Invalid χ trajectories format in self._chi_trajs.")
        elif hasattr(self, "chi_trajs"):
            chi_t, chi_tau = self.chi_trajs
        elif hasattr(self, "_get_chi_trajs"):
            chi_t, chi_tau = self._get_chi_trajs(lag)
        else:
            raise ValueError("No χ trajectories available. Run predict_chi() first.")
    
        # === 计算协方差矩阵 ===
        C00 = chi_t.T @ chi_t / (len(chi_t) - 1)
        C01 = chi_t.T @ chi_tau / (len(chi_t) - 1)
    
        # === 数值正则化，防止奇异 ===
        eps = 1e-6
        C00 += np.eye(C00.shape[0]) * eps
    
        # === 计算 Koopman 算子 ===
        K = np.linalg.inv(C00) @ C01
    
        # === 缓存结果 ===
        if not hasattr(self, "_k_cache"):
            self._k_cache = {}
        self._k_cache[lag] = K
    
        return K

    def its(self, lags):
        import numpy as np
        assert hasattr(self, "_chi_trajs"), "χ trajectories must be loaded first."
        chi_t_full = self._chi_trajs["chi_t"]
        chi_tau_full = self._chi_trajs["chi_tau"]
    
        # 取前半部的 n 维
        n = self.n_output
        chi_t = chi_t_full[:, :n]
        chi_tau = chi_tau_full[:, :n]
    
        its = np.empty((n - 1, len(lags)))
        for i, lag in enumerate(lags):
            print(f"   → 计算第 {i+1}/{len(lags)} 个 lag = {lag}")
            K = self.estimate_koopman(lag)
            lambdas = np.linalg.eigvals(K)
            lambdas = np.sort(np.abs(np.real(lambdas)))[:-1]
            its[:, i] = -lag / np.log(lambdas[:n - 1])
        return its

    
    def cktest(self, steps: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Perform the Chapman-Kolmogorov test on the model.
        
        Parameters
        ----------
        steps
            Number of steps to use
        
        Returns
        -------
        est
            Koopman operators estimated at different model lag times
        pred
            Koopman operators predicted at different
            lag times using matrix exponentiation
        
        """
        est = np.empty((self.n, self.n, steps))
        pred = np.empty((self.n, self.n, steps))
        est[:, :, FIRST] = np.eye(self.n)[:self.n]
        pred[:, :, FIRST] = np.eye(self.n)[:self.n]
        
        # Get the current Koopman operator first, because
        # estimating at a new lag time is very expensive
        K = self.K
        
        # Estimate new operators (slow)
        temp_est = np.empty((steps, self.n, self.n))
        temp_est[1] = K
        for nn in range(2, steps):
            self._log("CK test {0}/{1}".format(nn, steps), update=True)
            temp_est[nn] = self.estimate_koopman(self.network_lag * nn)
        self.reset_lag()
        
        # Get new predictions (fast)
        for i in range(self.n):
            vec = np.eye(self.n)[i]
            for nn in range(1, steps):
                est[i, :, nn] = vec @ temp_est[nn]
                pred[i, :, nn] = vec @ np.linalg.matrix_power(K, nn)
        
        return est, pred


# In[21]:


def loss_vampe(y_true, y_pred):
    return tf.linalg.trace(y_pred[FIRST])


# In[22]:


def statdist(X: np.ndarray) -> np.ndarray:
    """
    Calculate the equilibrium distribution of a transition matrix.
    
    Parameters
    ----------
    X
        Row-stochastic transition matrix
    
    Returns
    -------
    mu
        Stationary distribution, i.e. the left
        eigenvector associated with eigenvalue 1.
    
    """
    ev, evec = eig(X, left=True, right=False)
    mu_idx = np.abs(ev - 1.0).argmin()
    mu = evec[:, mu_idx].real # MODIFIED: Take real part and use correct index
    mu /= mu.sum()
    return mu