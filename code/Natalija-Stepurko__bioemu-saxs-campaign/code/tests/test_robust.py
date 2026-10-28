import numpy as np

from bsc import robust


def test_bootstrap_interval_covers_the_known_median():
    rng = np.random.default_rng(0)
    x = rng.lognormal(mean=np.log(1.2), sigma=0.3, size=200)
    lo, hi = robust.bootstrap_ci(x, np.median, n_boot=500, seed=1)
    assert lo < 1.2 < hi and hi - lo < 0.25
    # a share statistic: interval of a Bernoulli mean covers the true share and narrows with n
    b = rng.random(400) < 0.37
    lo, hi = robust.bootstrap_ci(b.astype(float), robust.share, n_boot=500, seed=2)
    assert lo < 0.37 < hi and hi - lo < 0.12
    assert robust.bootstrap_ci(np.array([np.nan]))[0] != robust.bootstrap_ci(np.array([np.nan]))[0]


def test_kmer_clusters_keep_near_duplicates_together():
    rng = np.random.default_rng(5)
    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    base = "".join(rng.choice(aa, 120))
    # 95% identity: six substitutions; a tagged construct; a truncation; and two unrelated chains
    mutant = list(base)
    for i in rng.choice(120, 6, replace=False):
        mutant[i] = "W" if mutant[i] != "W" else "A"
    mutant = "".join(mutant)
    tagged = "MGSSHHHHHHSSGLVPRGSH" + base
    truncated = base[20:100]
    other1, other2 = "".join(rng.choice(aa, 150)), "".join(rng.choice(aa, 90))
    seqs = [base, mutant, tagged, truncated, other1, other2]
    labs = ["base", "mutant", "tagged", "truncated", "o1", "o2"]
    pairs = robust.near_duplicate_pairs(labs, seqs)
    found = {frozenset((p["a"], p["b"])) for p in pairs}
    assert frozenset(("base", "mutant")) in found
    assert frozenset(("base", "tagged")) in found
    assert frozenset(("base", "truncated")) in found
    assert not any("o1" in p or "o2" in p for p in found)
    cl = robust.cluster_sequences(seqs, robust.GROUP_SIMILARITY)
    assert cl[0] == cl[1] == cl[2] == cl[3]
    assert len({cl[0], cl[4], cl[5]}) == 3
    # the near-duplicate threshold (Jaccard) joins the 95%-identity pair and the tagged construct
    cl = robust.cluster_sequences(seqs, robust.NEAR_DUPLICATE, measure=robust.jaccard)
    assert cl[0] == cl[1] == cl[2]
    assert len({cl[0], cl[4], cl[5]}) == 3
