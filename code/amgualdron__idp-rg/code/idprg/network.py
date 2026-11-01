"""The network: an MLP from the standardized sequence descriptors to the two-exponent distribution of Rg.

    ln S_m = b_S + NU_F ln N + o_1(x),    u = u_0 + o_2(x),    nu = 1/3 + (2/3) expit(u)
    ln A   = c . [1, nu, nu^2, 1/N, <lambda>, s+]

o = MLP(x) (GELU, last layer initialised at zero), x = the standardized inputs (descriptors.INPUTS). The amplitude law
c is learned with the MLP, started from least squares on the per-protein maximum-likelihood fits of the fit proteins;
b_S and u_0 are fixed at the mean of their law-constrained per-protein fits (data/md_fits.csv). Training minimises the
mean negative log-likelihood per frame of the MD frames (1,024 frames per protein redrawn each batch) with AdamW and a
cosine schedule; the last epoch is kept.
"""
import json
import pathlib

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import optax

from . import distribution as dist
from .data import NU_F, ROOT

MODEL = ROOT / "model"
LAW = ("1", "nu", "nu^2", "1/N", "mean_hydro", "scd_pos")


class Net(eqx.Module):
    mlp: eqx.nn.MLP
    law: jax.Array
    off: tuple = eqx.field(static=True)          # (b_S, u_0)

    def __init__(self, n_in, off, law, key, width=512, depth=2):
        mlp = eqx.nn.MLP(n_in, 2, width, depth, activation=jax.nn.gelu, key=key)
        last = mlp.layers[-1]
        self.mlp = eqx.tree_at(lambda m: (m.layers[-1].weight, m.layers[-1].bias), mlp,
                               (jnp.zeros_like(last.weight), jnp.zeros_like(last.bias)))
        self.law, self.off = jnp.asarray(law), tuple(float(v) for v in off)

    def params(self, z, raw):
        """z: standardized inputs (proteins, n_in); raw: (ln N, <lambda>, s+). Returns ln A, ln S_m, u."""
        o = jax.vmap(self.mlp)(z)
        lnN, lam, spos = raw[:, 0], raw[:, 1], raw[:, 2]
        lnSm = self.off[0] + NU_F * lnN + o[:, 0]
        u = self.off[1] + o[:, 1]
        nu = 1 / 3 + (2 / 3) * jax.nn.sigmoid(u)
        X = jnp.stack([jnp.ones_like(nu), nu, nu ** 2, jnp.exp(-lnN), lam, spos], -1)
        return X @ self.law, lnSm, u


def _raw(t):
    return t[["log_seq_length", "mean_hydro", "scd_pos"]].to_numpy(float)


def load(directory=MODEL):
    """(net, spec) of a saved network (default: the released model)."""
    directory = pathlib.Path(directory)
    spec = json.loads((directory / "network.json").read_text())
    net = Net(len(spec["inputs"]), (spec["b_S"], spec["u_0"]), np.zeros(len(LAW)), jax.random.PRNGKey(0),
              spec["width"], spec["depth"])
    return eqx.tree_deserialise_leaves(directory / "network.eqx", net), spec


def save(net, spec, directory):
    directory = pathlib.Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    eqx.tree_serialise_leaves(directory / "network.eqx", net)
    (directory / "network.json").write_text(json.dumps(spec, indent=1))


def parameters(net, spec, table):
    """(ln A, ln S_m, u) for the rows of a descriptor table (descriptors.table or data/proteins.csv)."""
    z = (table[spec["inputs"]].to_numpy(float) - np.array(spec["input_mean"])) / np.array(spec["input_sd"])
    return tuple(np.asarray(a) for a in net.params(jnp.asarray(z), jnp.asarray(_raw(table))))


def law_start(proteins, md_fits, fit_ids):
    """Least squares of the per-protein ln A on the law covariates, over the fit proteins."""
    t = md_fits.join(proteins.set_index("id")[["N", "mean_hydro", "scd_pos"]], on="id")
    m = t.id.isin(set(fit_ids)).to_numpy()
    nu = dist.nu_of(t.u.to_numpy())
    X = np.stack([np.ones_like(nu), nu, nu ** 2, 1 / t.N.to_numpy(float), t.mean_hydro.to_numpy(), t.scd_pos.to_numpy()], -1)
    c, *_ = np.linalg.lstsq(X[m], t.lnA.to_numpy()[m], rcond=None)
    return c


def train(proteins, md_fits, fit, F, law_split, inputs, width=512, depth=2, epochs=40, seed=42, log=print):
    """Train on the rows `fit` of the protein table. law_split names the md_fits columns (lnSm_<split>, u_<split>) of
    the law-constrained fits that set b_S and u_0. Returns (net, spec)."""
    t = proteins.join(md_fits.set_index("id")[[f"lnSm_{law_split}", f"u_{law_split}"]], on="id")
    X = t[list(inputs)].to_numpy(float)
    mu, sd = X[fit].mean(0), X[fit].std(0)
    Z, RAW = (X - mu) / sd, _raw(t)
    off = (float(np.mean(t[f"lnSm_{law_split}"].iloc[fit] - NU_F * t.log_seq_length.iloc[fit])),
           float(t[f"u_{law_split}"].iloc[fit].mean()))
    law0 = law_start(proteins, md_fits, t.id.iloc[fit])
    frames = np.asarray(F[t.frame_row.iloc[fit].to_numpy()], np.float32)

    net = Net(len(inputs), off, law0, jax.random.PRNGKey(seed), width, depth)
    nb = int(np.ceil(len(fit) / 64))
    opt = optax.chain(optax.clip_by_global_norm(1.0),
                      optax.adamw(optax.cosine_decay_schedule(1e-3, epochs * nb), weight_decay=1e-4))
    state = opt.init(eqx.filter(net, eqx.is_array))

    @eqx.filter_jit
    def step(n, st, z, raw, x):
        def loss(n):
            return -jnp.mean(dist.jax_logpdf(x, *n.params(z, raw)))
        l, g = eqx.filter_value_and_grad(loss)(n)
        upd, st = opt.update(g, st, eqx.filter(n, eqx.is_array))
        return eqx.apply_updates(n, upd), st, l

    rng = np.random.default_rng(seed)
    for ep in range(epochs):
        perm, ls = rng.permutation(len(fit)), []
        for b in range(0, len(fit), 64):
            idx = perm[b:b + 64]
            cols = rng.integers(0, frames.shape[1], size=(len(idx), 1024))
            x = np.log(np.take_along_axis(frames[idx], cols, axis=1).astype(float))
            net, state, l = step(net, state, jnp.asarray(Z[fit[idx]]), jnp.asarray(RAW[fit[idx]]), jnp.asarray(x))
            ls.append(float(l))
        log(f"epoch {ep + 1:3d}/{epochs}  training NLL {np.mean(ls):.5f}")
    spec = dict(inputs=list(inputs), input_mean=mu.tolist(), input_sd=sd.tolist(), law_covariates=list(LAW),
                law=np.asarray(net.law).tolist(), b_S=off[0], u_0=off[1], nu_F=NU_F, width=width, depth=depth,
                activation="gelu", epochs=epochs, seed=seed, n_fit=len(fit))
    return net, spec
