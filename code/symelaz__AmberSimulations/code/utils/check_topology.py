"""
check_topology.py <parm7> <MODEL> [WATER]

Checks that the intended force field really ended up in the topology
(tleap silently uses AmberTools files when a file with the same name exists;
this happened in the original pipeline_amber_tip4pd). Lines starting with FAIL
make build.sh stop.

MODEL: ff19SB_tip3p | a99SBdisp_tip4pd      WATER (a99SBdisp_tip4pd): authors | amber
"""
import sys
from collections import Counter
import parmed

parm = parmed.load_file(sys.argv[1])
MODEL = sys.argv[2]
WATER = sys.argv[3] if len(sys.argv) > 3 else ''
ok = True

def check(name, cond, detail):
    global ok
    print(f"{'PASS' if cond else 'FAIL'}  {name:40s} {detail}")
    ok &= bool(cond)

res = Counter(r.name for r in parm.residues)
print(f"MODEL {MODEL} {WATER} | atoms {len(parm.atoms)} | residues {len(parm.residues)} | box {parm.box[:3] if parm.box is not None else None}")
print("residue counts (> 20 or ions):", {k: v for k, v in res.items() if v > 20 or k in ('Na+', 'K+', 'Cl-')})

# --- water: (n_sites, H charge, O sigma A, O epsilon kcal/mol) ---
ref = {'tip3p':   (3, 0.417, 3.1507, 0.1520),
       'authors': (4, 0.59,  3.1650, 0.238764),
       'amber':   (4, 0.58,  3.1650, 0.238764)}[WATER]
water = [r for r in parm.residues if r.name in ('WAT', 'HOH')]
check("water present", len(water) > 0, f"{len(water)} waters")
if water:
    w = water[0]
    q = {a.name: round(a.charge, 4) for a in w.atoms}
    ow = [a for a in w.atoms if a.element == 8][0]
    check(f"water sites = {ref[0]}", len(w.atoms) == ref[0], str(q))
    check(f"water H charge = +{ref[1]}", abs(q['H1'] - ref[1]) < 1e-3, f"H1 {q['H1']}")
    check(f"water O LJ sigma {ref[2]}, eps {ref[3]}", abs(ow.sigma - ref[2]) < 2e-3 and abs(ow.epsilon - ref[3]) < 1e-4,
          f"sigma {ow.sigma:.4f} A, eps {ow.epsilon:.6f} kcal/mol  (TIP4P-Ew would be 3.1644 / 0.16275)")

types = Counter(a.type for a in parm.atoms)
if MODEL == 'ff19SB_tip3p':
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

q = sum(a.charge for a in parm.atoms)
check("net charge ~ 0", abs(q) < 1e-3, f"{q:.4f}")
print("RESULT:", "OK" if ok else "FAILED")
