"""
check_topology.py <parm7> <MODEL> [WATER] [disulfides.json]

Checks that the intended force field really ended up in the topology
(tleap silently uses AmberTools files when a file with the same name exists;
this happened in the original pipeline_amber_tip4pd). Lines starting with FAIL
make build.sh stop.

Always prints the water model's measured charges / LJ / O-EP bond length as
an "info" block, independent of pass/fail, so you can eyeball it against the
water model table in README.md and be sure which water you actually got.

MODEL: ff19SB_tip3p | ff19SB_tip4pew | a99SBdisp_tip4pd
WATER (a99SBdisp_tip4pd only): authors | amber
disulfides.json: written by utils/find_disulfides.py; if given and non-empty,
each recorded pair is checked for an actual SG-SG bond in the topology.
"""
import json
import sys
from collections import Counter

import parmed

parm = parmed.load_file(sys.argv[1])
MODEL = sys.argv[2]
WATER = sys.argv[3] if len(sys.argv) > 3 else ''
SS_JSON = sys.argv[4] if len(sys.argv) > 4 else None
ok = True


def check(name, cond, detail):
    global ok
    print(f"{'PASS' if cond else 'FAIL'}  {name:40s} {detail}")
    ok &= bool(cond)


res = Counter(r.name for r in parm.residues)
print(f"MODEL {MODEL} {WATER} | atoms {len(parm.atoms)} | residues {len(parm.residues)} | box {parm.box[:3] if parm.box is not None else None}")
print("residue counts (> 20 or ions):", {k: v for k, v in res.items() if v > 20 or k in ('Na+', 'K+', 'Cl-')})

# --- water: (n_sites, H charge, O sigma A, O epsilon kcal/mol, O-EP bond A or None) ---
# Values match the table in README.md ("Water models" section) - keep both in sync.
ref = {'tip3p':   (3, 0.417,   3.1507, 0.1520,    None),
       'tip4pew': (4, 0.52422, 3.1644, 0.16275,   0.1250),
       'authors': (4, 0.59,    3.1650, 0.238764,  0.1546),
       'amber':   (4, 0.58,    3.1650, 0.238764,  0.1546)}[WATER]
water = [r for r in parm.residues if r.name in ('WAT', 'HOH')]
check("water present", len(water) > 0, f"{len(water)} waters")
if water:
    w = water[0]
    q = {a.name: round(a.charge, 5) for a in w.atoms}
    ow = [a for a in w.atoms if a.element == 8][0]
    ep = next((a for a in w.atoms if a.name in ('EPW', 'EP')), None)
    d_om = None
    if ep is not None:
        b = next((b for b in ow.bonds if ep in (b.atom1, b.atom2)), None)
        d_om = round(b.type.req, 4) if b is not None and b.type is not None else None

    # Always-printed fingerprint: compare this block to the water model table
    # in README.md by eye, regardless of whether the PASS/FAIL checks below agree.
    print("--- water model fingerprint (compare to README.md) ---")
    print(f"  sites per molecule   : {len(w.atoms)}")
    print(f"  O  charge            : {q.get('O', q.get('OW', 'n/a'))}")
    print(f"  H  charge            : {q.get('H1', 'n/a')}")
    if ep is not None:
        print(f"  EP charge            : {q[ep.name]}")
        print(f"  O-EP bond length (A) : {d_om if d_om is not None else 'n/a'}")
    print(f"  O  sigma   (A)       : {ow.sigma:.4f}")
    print(f"  O  epsilon (kcal/mol): {ow.epsilon:.6f}")

    check(f"water sites = {ref[0]}", len(w.atoms) == ref[0], str(q))
    check(f"water H charge = +{ref[1]}", abs(q['H1'] - ref[1]) < 1e-3, f"H1 {q['H1']}")
    check(f"water O LJ sigma {ref[2]}, eps {ref[3]}", abs(ow.sigma - ref[2]) < 2e-3 and abs(ow.epsilon - ref[3]) < 1e-4,
          f"sigma {ow.sigma:.4f} A, eps {ow.epsilon:.6f} kcal/mol")
    if ref[4] is not None:
        check(f"water O-EP bond length {ref[4]}", d_om is not None and abs(d_om - ref[4]) < 2e-3,
              f"measured {d_om}")

types = Counter(a.type for a in parm.atoms)
if MODEL in ('ff19SB_tip3p', 'ff19SB_tip4pew'):
    # ff19SB uses residue-specific CMAP backbone corrections
    check("ff19SB CMAP terms present", len(parm.cmaps) > 0, f"{len(parm.cmaps)} CMAP terms")

if MODEL == 'a99SBdisp_tip4pd':
    pt = {t: types.get(t, 0) for t in ('C1', 'C8', 'C9', 'O3')}
    check("a99SB-disp atom types (C1/C8/C9/O3)", sum(pt.values()) > 0, f"{pt}  (0 everywhere = plain ff99SB-ILDN libraries)")
    o = next((a for a in parm.atoms if a.type == 'O'), None)
    h = next((a for a in parm.atoms if a.type == 'H'), None)
    if o is not None and h is not None:
        i, j = o.nb_idx - 1, h.nb_idx - 1
        acoef = parm.parm_data['LENNARD_JONES_ACOEF'][parm.parm_data['NONBONDED_PARM_INDEX'][parm.ptr('ntypes') * i + j] - 1]
        # LJEDIT "O H 2.2612 0.05742 0.0 0.05742": Rmin = 2.2612 + 0.0, eps = sqrt(0.05742*0.05742)
        # NOTE: this equals the normal combining rule for O/H, i.e. the LJEDIT changes nothing
        a_ref = 0.05742 * 2.2612 ** 12
        check("a99SB-disp LJEDIT O-H pair", abs(acoef / a_ref - 1) < 1e-3, f"A = {acoef:.1f} (expected {a_ref:.1f})")
    for name, r in (('K+', (1.7638, 0.087)), ('Cl-', (2.27, 0.15))):
        ion = next((a for a in parm.atoms if a.residue.name == name), None)
        if ion is None:
            print(f"info  no {name} in system")
            continue
        check(f"{name} CHARMM22 LJ", abs(ion.rmin - r[0]) < 1e-3 and abs(ion.epsilon - r[1]) < 1e-4,
              f"Rmin/2 {ion.rmin:.4f}, eps {ion.epsilon:.4f} (ref {r})")
    n_extra = sum(1 for d in parm.dihedrals if d.type is not None and abs(d.type.phase - 90.0) < 1e-6
                  and abs(d.type.phi_k - 0.40) < 1e-6 and d.type.per == 2)
    check("add8 backbone dihedrals present", n_extra > 0, f"{n_extra} phi terms (k=0.40, n=2, 90 deg)")

# --- disulfides: each pair find_disulfides.py recorded must be an actual SG-SG
# bond in the final topology, both residues CYX, and no stray HG re-added. ---
if SS_JSON:
    try:
        pairs = json.load(open(SS_JSON))
    except FileNotFoundError:
        pairs = []
    if pairs:
        print(f"--- disulfides (from {SS_JSON}) ---")
        for p in pairs:
            i, j = p['leap_index']
            ri, rj = parm.residues[i - 1], parm.residues[j - 1]
            sgi = next((a for a in ri.atoms if a.name == 'SG'), None)
            sgj = next((a for a in rj.atoms if a.name == 'SG'), None)
            bonded = sgi is not None and sgj is not None and any(
                sgj in (b.atom1, b.atom2) for b in sgi.bonds)
            no_hg = not any(a.name == 'HG' for a in ri.atoms) and not any(a.name == 'HG' for a in rj.atoms)
            is_cyx = ri.name == 'CYX' and rj.name == 'CYX'
            label = f"{p['chain_resnum'][0]}-{p['chain_resnum'][1]} (input SG-SG {p['distance_A']} A)"
            check(f"disulfide {label}", bonded and no_hg and is_cyx,
                  f"SG-SG bond={bonded}, both CYX={is_cyx}, no stray HG={no_hg}")
    else:
        print("info  no disulfides were detected in the input (0 CYS pairs within cutoff)")

# Every CYX must actually be bonded to another SG: a CYX without its partner
# has no HG either, i.e. a sulfur radical (happens when the input already
# names a residue CYX but its partner is out of range or missing).
lonely = []
for r in parm.residues:
    if r.name != 'CYX':
        continue
    sg = next((a for a in r.atoms if a.name == 'SG'), None)
    if sg is None or not any(o.name == 'SG' and o.residue is not r
                             for b in sg.bonds for o in (b.atom1, b.atom2)):
        lonely.append(f"{r.name}{r.idx + 1}")
check("every CYX has an SG-SG bond", not lonely, f"unbonded: {lonely}" if lonely else f"{sum(r.name == 'CYX' for r in parm.residues)} CYX")

# Peptide bonds: tleap bonds C(i)-N(i+1) across any break it is not told about
# with a TER, printing only a warning. Check the built coordinates.
try:
    import numpy as np
    xyz = parmed.load_file(sys.argv[1], xyz=sys.argv[1].rsplit('.', 1)[0] + '.rst7').coordinates
    long_cn = []
    for b in parm.bonds:
        a1, a2 = b.atom1, b.atom2
        if a1.residue is not a2.residue and {a1.name, a2.name} == {'C', 'N'}:
            d = float(np.linalg.norm(xyz[a1.idx] - xyz[a2.idx]))
            if d > 2.0:
                long_cn.append(f"{a1.residue.name}{a1.residue.idx + 1}-{a2.residue.name}{a2.residue.idx + 1}:{d:.2f}A")
    check("no stretched peptide bonds (C-N <= 2 A)", not long_cn,
          " ".join(long_cn) if long_cn else "all inter-residue C-N bonds <= 2 A")
    # Any other bond (ligand included; a disulfide is ~2.05 A): a wrong
    # template/connectivity shows up as a bond between atoms far apart.
    long_any = []
    for b in parm.bonds:
        d = float(np.linalg.norm(xyz[b.atom1.idx] - xyz[b.atom2.idx]))
        if d > 2.5:
            long_any.append(f"{b.atom1.residue.name}{b.atom1.residue.idx + 1}:{b.atom1.name}-"
                            f"{b.atom2.residue.name}{b.atom2.residue.idx + 1}:{b.atom2.name}:{d:.2f}A")
    check("no stretched bonds anywhere (<= 2.5 A)", not long_any,
          " ".join(long_any[:10]) if long_any else "all bonds <= 2.5 A")
except FileNotFoundError:
    print("info  no .rst7 next to the topology, peptide bond lengths not checked")

q = sum(a.charge for a in parm.atoms)
check("net charge ~ 0", abs(q) < 1e-3, f"{q:.4f}")
print("RESULT:", "OK" if ok else "FAILED")
