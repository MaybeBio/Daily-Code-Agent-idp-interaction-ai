"""
design_sequences.py — the GAIVI design batch for the RECOMB hub-propensity paper: designed sequences plus the
experimental PED constructs, simulated in CALVADOS-2, Mpipi-Recharged and MPIPI-GG with the BENDER protocol.
These are prospective tests: predictions are frozen before any of the labels are computed.

Writes next to this file:
  hub_design_300K.fasta   >ID|category|key=value;...   (the ID is the run name in all three submit scripts;
                          "300K" in the file name sets the temperature for submit_calvados.sh / submit_openmpipi.sh)
  hub_design_table.csv    id, category, n, parent, host, guest, positions, notes, sequence

Categories (hosts and designs are 100 residues):
  ped          the 14 distinct sequences of the 18 PED ensembles, exact constructs (tags included). pSic1 is the
               phosphomimetic (pS/pT -> E), since none of the three models has phosphoresidues.
  variant      for the 13 unmodified PED constructs: aromatics -> S, R -> K, K -> R, D/E -> N/Q, two shuffles
               (composition kept, order destroyed). Only variants that change the sequence are made.
  homopolymer  the 20 amino acids x 100: purely positional profiles, no residue-type contrast is possible
  host         GS = (GS)50; AVG = mean BENDER composition (aromatics included), random order; CHG = E/K-rich
               polyampholyte, well mixed; PLD = prion-like S/G/Q/N/T, aromatics removed
  guest        every amino acid at the centre (index 50) of each host; Y, W, R, K, E, P also at indices 10 and 25
               of GS (position x residue)
  valence      Y, W, R, K, E in GS: 2 / 4 / 8 / 12 copies evenly spaced, or 2 / 4 / 8 in one central block
  pair         one Y plus one partner (R, K, E, Y) in GS at sequence separation 1, 2, 4, 8, 16, 32
  kappa        E/K 50:50 with charge patterning kappa 0 ... 1 (Das & Pappu 2013 definition, blobs 5 and 6)
  replica      12 BENDER sequences already simulated in all three force fields (2 per group, N <= 150, from the
               908 paired set): a second, independent run gives the between-run noise ceiling of the labels
  (appended)   8 more 'ped' constructs never used before (fresh external test; ../ped_fresh_set.csv) and 4 'pre'
               constructs: the disordered proteins with PRE data in PeptoneDB-Integrative (OPN, tau 2N4R, ER NTD,
               ANAC046), whose simulations will be reweighted against that PRE
Deterministic (fixed seeds). Standard library only.  python3 design_sequences.py
"""
import csv
import os
import random
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
sys.path.insert(0, EV)
AAS = 'ACDEFGHIKLMNPQRSTVWY'
N = 100
CENTRE = 50
BENDER_FASTA = os.path.expanduser('~/Downloads/Scripts/OpenMPIPI/bender_10282.fasta')
PEPTONE = os.path.expanduser('~/Downloads/PeptoneDBs/PeptoneDB-Integrative')      # PeptoneBench data (Zenodo 17306061)
PED_NAMES = {'PED00024': 'aSyn', 'PED00003': 'bSyn', 'PED00006': 'abSyn', 'PED00022': 'drkN', 'PED00023': 'Sic1',
             'PED00161': 'pSic1E', 'PED00225': 'eIF4G1', 'PED00498': 'V2R', 'PED00499': 'GHSR', 'PED00500': 'B2AR',
             'PED00511': 'PDEg', 'PED00013': 'PaaA2', 'PED00004': 'p27KID', 'PED00494': 'FUSNTD'}
rows = []


def add(uid, cat, seq, parent='', host='', guest='', positions='', notes=''):
    assert set(seq) <= set(AAS), (uid, set(seq) - set(AAS))
    rows.append({'id': uid, 'category': cat, 'n': len(seq), 'parent': parent, 'host': host, 'guest': guest,
                 'positions': positions, 'notes': notes, 'sequence': seq})


def from_counts(counts, seed):
    s = [a for a, k in counts.items() for _ in range(k)]
    assert len(s) == N, len(s)
    random.Random(seed).shuffle(s)
    return ''.join(s)


def round_counts(freq, total=N):
    """Largest-remainder rounding of total * freq to integers that sum to total."""
    raw = {a: total * f / sum(freq.values()) for a, f in freq.items()}
    k = {a: int(v) for a, v in raw.items()}
    for a in sorted(raw, key=lambda a: -(raw[a] - k[a]))[:total - sum(k.values())]:
        k[a] += 1
    return k


def mutate(seq, pos, aa):
    s = list(seq)
    for p in pos:
        s[p] = aa
    return ''.join(s)


# ── charge patterning ────────────────────────────────────────────────────────
def _delta(s, g):
    n = len(s)
    fp, fm = sum(c in 'KR' for c in s) / n, sum(c in 'DE' for c in s) / n
    sig = (fp - fm) ** 2 / (fp + fm) if fp + fm else 0.0
    tot = 0.0
    for i in range(n - g + 1):
        b = s[i:i + g]
        p, m = sum(c in 'KR' for c in b) / g, sum(c in 'DE' for c in b) / g
        tot += (((p - m) ** 2 / (p + m) if p + m else 0.0) - sig) ** 2
    return tot / (n - g + 1)


def kappa(s):
    """kappa of a sequence whose charged residues are K/R and D/E (for pure E/K the maximum is the two-block chain)."""
    seg = ''.join(sorted(s, key=lambda c: 0 if c in 'KR' else 2 if c in 'DE' else 1))
    return sum(_delta(s, g) / _delta(seg, g) for g in (5, 6)) / 2


def kappa_design(target, seed, steps=40000):
    rng = random.Random(seed)
    s = list('EK' * (N // 2))
    rng.shuffle(s)
    cur = abs(kappa(''.join(s)) - target)
    for _ in range(steps):
        i, j = rng.randrange(N), rng.randrange(N)
        if s[i] == s[j]:
            continue
        s[i], s[j] = s[j], s[i]
        new = abs(kappa(''.join(s)) - target)
        if new <= cur or rng.random() < 0.01:
            cur = new
        else:
            s[i], s[j] = s[j], s[i]
        if cur < 0.005:
            break
    return ''.join(s)


def main():
    import ped_null as pn
    # ped + variants
    seen = {}
    for e, f, prot, _ in pn.ENSEMBLES:
        seq, _ = pn.load(f)
        if seq in seen:
            continue
        seen[seq] = e
        name = PED_NAMES[e]
        if e == 'PED00161':
            add(f'PED_{name}', 'ped', seq.replace('x', 'E').replace('y', 'E').replace('z', 'E'), parent=e,
                notes=f'{prot}; phosphomimetic: {sum(c in "xyz" for c in seq)} phosphoresidues -> E')
            continue
        add(f'PED_{name}', 'ped', seq, parent=e, notes=prot)
        for tag, src, dst, what in (('aro2S', 'FWY', 'SSS', 'aromatics -> S'), ('R2K', 'R', 'K', 'R -> K'),
                                    ('K2R', 'K', 'R', 'K -> R'), ('DE2NQ', 'DE', 'NQ', 'D -> N, E -> Q')):
            v = seq.translate(str.maketrans(src, dst))
            if v != seq:
                pos = [i for i, (a, b) in enumerate(zip(seq, v)) if a != b]
                add(f'VAR_{name}_{tag}', 'variant', v, parent=f'PED_{name}', positions=' '.join(map(str, pos)),
                    notes=f'{what} ({len(pos)} sites)')
        for k in (1, 2):
            s = list(seq)
            random.Random(zlib.crc32(f'{name}{k}'.encode())).shuffle(s)
            add(f'VAR_{name}_shuf{k}', 'variant', ''.join(s), parent=f'PED_{name}', notes='shuffle (same composition)')

    # homopolymers
    for a in AAS:
        add(f'HOMO_{a}100', 'homopolymer', a * N, guest=a)

    # hosts
    comp = {a: 0 for a in AAS}
    for line in open(BENDER_FASTA):
        if not line.startswith('>'):
            for c in line.strip():
                if c in comp:
                    comp[c] += 1
    hosts = {
        'GS': 'GS' * (N // 2),
        'AVG': from_counts(round_counts(comp), 101),
        'CHG': from_counts({'E': 25, 'K': 25, 'S': 15, 'G': 10, 'P': 10, 'A': 10, 'Q': 5}, 102),
        'PLD': from_counts({'S': 30, 'G': 25, 'Q': 20, 'T': 8, 'N': 6, 'A': 5, 'P': 4, 'D': 2}, 103),
    }
    notes = {'GS': '(GS)50', 'AVG': 'mean BENDER composition, random order', 'CHG': 'E/K-rich polyampholyte, random order',
             'PLD': 'prion-like S/G/Q/N/T, no aromatics, random order'}
    for h, s in hosts.items():
        add(f'HOST_{h}', 'host', s, host=h, notes=notes[h])
        for a in AAS:
            if s[CENTRE] != a:
                add(f'GUEST_{h}_{a}{CENTRE}', 'guest', mutate(s, [CENTRE], a), parent=f'HOST_{h}', host=h, guest=a,
                    positions=str(CENTRE), notes=f'{s[CENTRE]}{CENTRE}{a}')
    for a in 'YWRKEP':
        for p in (10, 25):
            add(f'GUEST_GS_{a}{p}', 'guest', mutate(hosts['GS'], [p], a), parent='HOST_GS', host='GS', guest=a,
                positions=str(p), notes=f'{hosts["GS"][p]}{p}{a}')

    # valence and patterning of one residue type
    for a in 'YWRKE':
        for k in (2, 4, 8, 12):
            pos = [int((j + 0.5) * N / k) for j in range(k)]
            add(f'VAL_{a}_even{k}', 'valence', mutate(hosts['GS'], pos, a), parent='HOST_GS', host='GS', guest=a,
                positions=' '.join(map(str, pos)), notes=f'{k} x {a}, evenly spaced')
        for k in (2, 4, 8):
            pos = list(range(CENTRE - k // 2, CENTRE - k // 2 + k))
            add(f'VAL_{a}_block{k}', 'valence', mutate(hosts['GS'], pos, a), parent='HOST_GS', host='GS', guest=a,
                positions=' '.join(map(str, pos)), notes=f'{k} x {a}, one central block')

    # sticker pairs
    for b in 'RKEY':
        for d in (1, 2, 4, 8, 16, 32):
            i = CENTRE - d // 2
            add(f'PAIR_Y{b}_d{d}', 'pair', mutate(mutate(hosts['GS'], [i], 'Y'), [i + d], b), parent='HOST_GS',
                host='GS', guest=f'Y+{b}', positions=f'{i} {i + d}', notes=f'Y{i} and {b}{i + d}, separation {d}')

    # charge patterning
    for t in (0.0, 0.1, 0.2, 0.3, 0.4, 0.55, 0.7, 1.0):
        if t == 0.0:
            s = 'EK' * (N // 2)
        elif t == 1.0:
            s = 'E' * (N // 2) + 'K' * (N // 2)
        else:
            s = kappa_design(t, int(t * 1000))
        add(f'KAPPA_EK_{int(round(t * 100)):03d}', 'kappa', s, notes=f'E/K 50:50, kappa {kappa(s):.3f} (target {t})')

    # replicas of BENDER sequences simulated in all three force fields
    fa, cur = {}, None
    for line in open(BENDER_FASTA):
        line = line.strip()
        if line.startswith('>'):
            cur = line[1:].split('|')[0]
            fa[cur] = ''
        elif cur:
            fa[cur] += line
    paired = list(csv.DictReader(open(os.path.join(EV, 'hub_frames_mpipi', 'ids.csv'))))
    rech = {r['UniProt_ID'] for r in csv.DictReader(open(os.path.join(EV, 'hub_frames_recharged', 'per_protein.csv')))}
    by = {}
    for r in sorted(paired, key=lambda r: r['UniProt_ID']):
        if r['UniProt_ID'] in fa and len(fa[r['UniProt_ID']]) <= 150:
            by.setdefault(r['group'], []).append(r['UniProt_ID'])
    rng = random.Random(7)
    for g in sorted(by):
        ids = by[g][:]
        rng.shuffle(ids)
        ids.sort(key=lambda u: u not in rech)                # prefer IDs that already have Recharged labels
        for u in ids[:2]:
            add(f'REP_{u}', 'replica', fa[u], parent=u, notes=f'{g}; replica of the BENDER run'
                + ('; Recharged pilot label exists' if u in rech else ''))

    # added 2026-10-02 (appended, so the first 280 rows match predictions_frozen_2026-10-02.csv):
    # 8 fresh PED constructs never used before (disorder-checked; predictions frozen in *_ped_fresh.csv) ...
    fresh = {'PED00501': 'V2Rpm', 'PED00502': 'GHSRpm', 'PED00503': 'B2ARpm', 'PED00154': 'MKK4', 'PED00017': 'TauF',
             'PED00533': 'Ab40', 'PED00216': 'CBPID4', 'PED00536': 'ACTR'}
    for r in csv.DictReader(open(os.path.join(EV, 'ped_fresh_set.csv'))):
        add(f"PED_{fresh[r['id']]}", 'ped', r['sequence'], parent=r['id'],
            notes=f"{r['protein']}; {r['restraints']}; fresh external test")
    # ... and the 4 disordered proteins with PRE data in PeptoneDB-Integrative that we have no ensemble for
    # (exact constructs from seq.fasta; their PRE is used to reweight these simulations)
    for p, name in (('OPN', 'OPN'), ('Tau', 'Tau441'), ('ERNTD', 'ERNTD'), ('ANAC046', 'ANAC046')):
        d = os.path.join(PEPTONE, p)
        seq = ''.join(line.strip() for line in open(os.path.join(d, 'seq.fasta')) if not line.startswith('>'))
        n_pre = sum(f.startswith('PRE') for f in os.listdir(d))
        add(f'PRE_{name}', 'pre', seq, parent=f'PeptoneDB-Integrative/{p}', notes=f'{n_pre} PRE spin-label sites')

    # optional tier 4 (appended): two more independent runs of every host and guest sequence, because single-run noise
    # in one residue's occupancy (SE >= ~0.009) is comparable to the predicted guest effects (median 0.011)
    for r in [r for r in rows if r['category'] in ('host', 'guest')]:
        for k in (2, 3):
            add(f"{r['id']}_r{k}", 'gh_replica', r['sequence'], parent=r['id'], host=r['host'], guest=r['guest'],
                positions=r['positions'], notes=f'independent run {k} of {r["id"]}')

    ids = [r['id'] for r in rows]
    assert len(ids) == len(set(ids)), 'duplicate IDs'
    seqs = {}
    for r in rows:
        if r['sequence'] in seqs and r['category'] != 'replica':
            print(f"note: {r['id']} has the same sequence as {seqs[r['sequence']]}")
        seqs.setdefault(r['sequence'], r['id'])
    with open(os.path.join(HERE, 'hub_design_300K.fasta'), 'w') as fh:
        for r in rows:
            meta = ';'.join(f'{k}={r[k]}' for k in ('parent', 'host', 'guest') if r[k])
            fh.write(f">{r['id']}|{r['category']}|N={r['n']}{';' + meta if meta else ''}\n{r['sequence']}\n")
    with open(os.path.join(HERE, 'hub_design_table.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['id', 'category', 'n', 'parent', 'host', 'guest', 'positions', 'notes',
                                           'sequence'])
        w.writeheader()
        w.writerows(rows)
    cats = {}
    for r in rows:
        cats[r['category']] = cats.get(r['category'], 0) + 1
    print(f'{len(rows)} sequences, {sum(r["n"] for r in rows)} residues:', cats)


if __name__ == '__main__':
    main()
