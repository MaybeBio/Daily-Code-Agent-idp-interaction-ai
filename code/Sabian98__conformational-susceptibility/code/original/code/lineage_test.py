"""
Does susceptibility at matched compaction differ across lineages, and is that more than composition?

Self-contained (no trajectory or NAS module access): reads two tables from --inputs
  per_protein.csv   energy_mechanism_split output: chi_C (% of <Rg> per K, mean over 278-320 K), nu_C, length_C, rg300_C
  BENDER_BIO.csv    organism, taxonomy_id, kingdom, sequence
CALVADOS-2, all chains with five-temperature chi.

Leftover L: chi residualised on a cubic spline of nu (5 knots) and log N, as in Test 2b; in SD units.
Composition-adjusted leftover L_adj: L residualised on rank-normalised pi fraction (F, W, Y, R) and D + E + K fraction.
Residue-age score: mean rank of a protein's residues in the proposed order of entry into the genetic code
(Kulkarni et al., Chem. Rev. 2022, section 5.1): G/A 1, V/D 2, P 3, S 4, E/L 5, T 6, R 7, N 8, K 9, Q 10, I 11, C 12,
H 13, F 14, M 15, Y 16, W 17 (higher = later residues).

Proteins are clustered in organisms (rice alone is 16% of chains), so lineage comparisons use organism means
(each organism weighs the same): kingdom = mean of its organisms' means, 95% CI by bootstrap over organisms,
kingdom effect by permutation of kingdom labels over organisms (one-way F). A mixed model (L ~ kingdom, random
intercept per organism) is reported alongside. Cellularity at kingdom level: unicellular = Bacteria + Protists,
multicellular = Mammals + Plants; Fungi (yeasts and filamentous) and Viruses kept separate.
Output: --out (default Traj Level Analysis/Analysis_FDT/lineage_test/)
Usage: python3 lineage_test.py --inputs DIR [--out DIR] [--model C|M] [--common]
(--model M or --common write to a sub-folder of --out)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

OUT_DEFAULT = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/lineage_test'
AGE = {**dict.fromkeys('GA', 1), **dict.fromkeys('VD', 2), 'P': 3, 'S': 4, **dict.fromkeys('EL', 5), 'T': 6, 'R': 7,
       'N': 8, 'K': 9, 'Q': 10, 'I': 11, 'C': 12, 'H': 13, 'F': 14, 'M': 15, 'Y': 16, 'W': 17}
CELL = {'Bacteria': 'unicellular', 'Protists': 'unicellular', 'Mammals': 'multicellular', 'Plants': 'multicellular',
        'Fungi': 'fungi (mixed)', 'Viruses': 'viruses'}
ORDER = ['Bacteria', 'Protists', 'Fungi', 'Plants', 'Mammals', 'Viruses']
N_BOOT, N_PERM = 2000, 5000


def spline(x, knots=5):
    q = np.quantile(x, np.linspace(0, 1, knots + 2)[1:-1])
    return np.column_stack([x, x ** 2, x ** 3] + [np.clip(x - t, 0, None) ** 3 for t in q])


def residualise(y, X):
    X1 = np.column_stack([np.ones(len(y)), X])
    return y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]


def rank_normal(x):
    return stats.norm.ppf((stats.rankdata(x) - 0.5) / len(x))


def frac(seq, residues):
    return seq.map(lambda q: sum(q.count(a) for a in residues) / len(q))


def load(inputs, model='C', common=False):
    cols = [f'{c}_{m}' for m in 'CM' for c in ('chi', 'nu', 'length')]
    per = pd.read_csv(os.path.join(inputs, 'per_protein.csv'), index_col=0, usecols=['UniProt_ID'] + cols)
    if common:
        per = per.dropna()
    per = per[[f'chi_{model}', f'nu_{model}', f'length_{model}']].set_axis(['chi_C', 'nu_C', 'length_C'], axis=1)
    per.index.name = 'UniProt_ID'
    bio = pd.read_csv(os.path.join(inputs, 'BENDER_BIO.csv'), usecols=['UniProt_ID', 'organism', 'taxonomy_id', 'kingdom', 'sequence'],
                      low_memory=False).set_index('UniProt_ID')
    d = per.dropna().join(bio, how='inner')
    d['pi'] = frac(d.sequence, 'FWYR')
    d['dek'] = frac(d.sequence, 'DEK')
    d['age'] = d.sequence.map(lambda q: np.mean([AGE[c] for c in q if c in AGE]))
    logN = np.log(d.length_C.to_numpy())
    L = residualise(d.chi_C.to_numpy(), np.column_stack([spline(d.nu_C.to_numpy()), logN]))
    d['L'] = L / L.std()
    La = residualise(d.L.to_numpy(), np.column_stack([rank_normal(d.pi.to_numpy()), rank_normal(d.dek.to_numpy())]))
    d['L_adj'] = La / d.L.std()
    d['cell'] = d.kingdom.map(CELL)
    return d


def org_means(d, cols):
    return d.groupby('taxonomy_id').agg(kingdom=('kingdom', 'first'), cell=('cell', 'first'), n=('L', 'size'),
                                        **{c: (c, 'mean') for c in cols})


def boot_ci(x, rng):
    b = rng.choice(x, (N_BOOT, len(x))).mean(1)
    return np.percentile(b, [2.5, 97.5])


def perm_F(values, groups, rng):
    def F(g):
        return stats.f_oneway(*[values[g == k] for k in np.unique(g)]).statistic
    f0 = F(groups)
    null = np.array([F(rng.permutation(groups)) for _ in range(N_PERM)])
    return f0, (np.sum(null >= f0) + 1) / (N_PERM + 1)


def perm_diff(a, b, rng):
    x = np.concatenate([a, b])
    d0 = a.mean() - b.mean()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        p = rng.permutation(x)
        null[i] = p[:len(a)].mean() - p[len(a):].mean()
    return d0, (np.sum(np.abs(null) >= abs(d0)) + 1) / (N_PERM + 1)


def main(inputs, out, model='C', common=False):
    out = out if (model == 'C' and not common) else os.path.join(out, f"{'mpipi' if model == 'M' else 'calvados'}{'_common' if common else ''}")
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(0)
    d = load(inputs, model, common)
    print(f"model {'Mpipi-GG' if model == 'M' else 'CALVADOS-2'}{', common set' if common else ''}")
    cols = ['L', 'L_adj', 'chi_C', 'nu_C', 'pi', 'dek', 'age']
    o = org_means(d, cols)
    print(f'{len(d)} chains, {len(o)} organisms; chains per kingdom {d.kingdom.value_counts().to_dict()}')
    print(f'largest organisms: {d.organism.value_counts().head(3).to_dict()}')

    # kingdoms: organism-weighted means with bootstrap CIs, raw and composition-adjusted
    rows = []
    for k in ORDER:
        x = o[o.kingdom == k]
        row = {'kingdom': k, 'cell': CELL[k], 'chains': int(d.kingdom.eq(k).sum()), 'organisms': len(x)}
        for c in cols:
            row[f'{c}_mean'] = x[c].mean()
        for c in ('L', 'L_adj'):
            row[f'{c}_lo'], row[f'{c}_hi'] = boot_ci(x[c].to_numpy(), rng)
        row['L_chain_median'] = d.loc[d.kingdom == k, 'L'].median()
        rows.append(row)
    K = pd.DataFrame(rows).set_index('kingdom')
    tests = []
    for c in ('L', 'L_adj', 'chi_C'):
        f, p = perm_F(o[c].to_numpy(), o.kingdom.to_numpy(), rng)
        tests.append({'test': f'kingdom effect on {c} (organism means, permutation F)', 'stat': f, 'p': p})
    cellular = o[o.cell.isin(['unicellular', 'multicellular'])]
    for c in ('L', 'L_adj', 'age'):
        diff, p = perm_diff(cellular.loc[cellular.cell == 'multicellular', c].to_numpy(),
                            cellular.loc[cellular.cell == 'unicellular', c].to_numpy(), rng)
        tests.append({'test': f'multicellular minus unicellular, {c} (organism means, permutation)', 'stat': diff, 'p': p})
    # residue age vs leftover, chain level and organism level; partial on composition
    for lvl, x in (('chains', d), ('organisms', o)):
        r, p = stats.spearmanr(x.age, x.L)
        tests.append({'test': f'Spearman residue-age score vs L ({lvl})', 'stat': r, 'p': p})
    ra = residualise(rank_normal(d.age.to_numpy()), np.column_stack([rank_normal(d.pi.to_numpy()), rank_normal(d.dek.to_numpy())]))
    r, p = stats.pearsonr(ra, d.L_adj)
    tests.append({'test': 'residue-age score vs L, both adjusted for pi and D+E+K (chains, rank-based partial)', 'stat': r, 'p': p})
    # mixed model: L ~ kingdom + (1 | organism)
    try:
        import statsmodels.formula.api as smf
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for c in ('L', 'L_adj'):
                m = smf.mixedlm(f'{c} ~ C(kingdom, Treatment("Bacteria"))', d.reset_index(), groups='taxonomy_id').fit(reml=True)
                for name, val in m.params.items():
                    if 'kingdom' in name:
                        tests.append({'test': f'mixed model {c}: {name.split("[T.")[1].rstrip("]")} minus Bacteria', 'stat': val,
                                      'p': m.pvalues[name]})
                tests.append({'test': f'mixed model {c}: organism variance share',
                              'stat': float(m.cov_re.iloc[0, 0] / (m.cov_re.iloc[0, 0] + m.scale)), 'p': np.nan})
    except Exception as e:  # noqa: BLE001
        tests.append({'test': f'mixed model failed: {e!r}', 'stat': np.nan, 'p': np.nan})
    T = pd.DataFrame(tests)
    d[['organism', 'taxonomy_id', 'kingdom', 'cell', 'chi_C', 'nu_C', 'length_C', 'pi', 'dek', 'age', 'L', 'L_adj']].to_csv(os.path.join(out, 'per_protein.csv'))
    o.to_csv(os.path.join(out, 'per_organism.csv'))
    K.to_csv(os.path.join(out, 'kingdoms.csv'))
    T.to_csv(os.path.join(out, 'tests.csv'), index=False)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    print(K.round(3).to_string())
    print(T.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    print(f'saved {out}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', required=True)
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--model', choices=['C', 'M'], default='C', help='C = CALVADOS-2, M = Mpipi-GG (columns chi/nu/length_<model>)')
    ap.add_argument('--common', action='store_true', help='only proteins simulated in both models')
    a = ap.parse_args()
    main(a.inputs, a.out, a.model, a.common)
