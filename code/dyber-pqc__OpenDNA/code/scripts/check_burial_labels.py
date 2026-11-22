"""Are the burial labels actually right?

    python scripts/check_burial_labels.py datasets/out/burial/data/train.parquet


If the sequence-to-coordinate alignment were off, the numbers would still look
plausible — a shifted label set has the same distribution. What a shift
destroys is the *correlation with chemistry*: buried positions must be
enriched in hydrophobic residues and depleted in charged ones. Checking that
tests the alignment end to end, which no unit test on the parser can.
"""

import pyarrow.parquet as pq

BURIAL = (
    r"H:\opendna\.claude\worktrees\relaxed-feistel-a92667\datasets\out"
    r"\burial\data\train.parquet"
)

# Kyte-Doolittle sign: these are the unambiguous cases either way.
HYDROPHOBIC = set("VILFMCWA")
CHARGED = set("KEDR")

t = pq.read_table(BURIAL).to_pydict()
print(f"domains: {len(t['domain_id'])}")

mean_rsa, counts = {}, {}
for seq, rsa in zip(t["sequence"], t["rsa"]):
    for aa, value in zip(seq, rsa):
        mean_rsa[aa] = mean_rsa.get(aa, 0.0) + value
        counts[aa] = counts.get(aa, 0) + 1

rows = sorted(((mean_rsa[a] / counts[a], a, counts[a]) for a in counts))
print("\nmean relative accessibility per residue, most buried first:")
for value, aa, n in rows:
    tag = "hydrophobic" if aa in HYDROPHOBIC else ("charged" if aa in CHARGED else "")
    print(f"  {aa}  {value:.3f}  n={n:<6} {tag}")

hyd = [v for v, a, _ in rows if a in HYDROPHOBIC]
chg = [v for v, a, _ in rows if a in CHARGED]
mh, mc = sum(hyd) / len(hyd), sum(chg) / len(chg)
print(f"\nmean rsa, hydrophobic residues: {mh:.3f}")
print(f"mean rsa, charged residues:     {mc:.3f}")
print(f"difference:                     {mc - mh:+.3f}")

top5 = [a for _, a, _ in rows[:5]]
bottom5 = [a for _, a, _ in rows[-5:]]
print(f"\nfive most buried:  {top5}")
print(f"five most exposed: {bottom5}")

assert mc > mh, "charged residues must be more exposed than hydrophobic ones"
assert sum(a in HYDROPHOBIC for a in top5) >= 3, f"expected mostly hydrophobic buried, got {top5}"
assert sum(a in CHARGED for a in bottom5) >= 2, f"expected charged exposed, got {bottom5}"
print("\nchemistry agrees with the labels: the alignment holds end to end")
