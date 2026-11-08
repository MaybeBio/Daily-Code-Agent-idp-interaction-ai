"""test_helpers.py: the released helpers/ must give exactly what the BENDER-MT originals give.
Runs one process per module set, on a real CALVADOS-2 trajectory and on a sequence with every residue type, and compares every output.
  BENDER_REANALYSIS_DIR=... ~/.venvs/finches/bin/python tools/test_helpers.py [UniProt_ID ...]
"""
import os
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
MT = "/Volumes/User Homes/Nu Project/BENDER-MT/code"
SEQ = "ACDEFGHIKLMNPQRSTVWYGGSSKKEEDDRRFFYYWW"
CHILD = r"""
import pickle, sys, os
sys.path.insert(0, sys.argv[1])
import numpy as np
import salt_ph_phospho as sp
out = {'file': sp.__file__, 'CAL': sp.CAL, 'BETA_C': sp.BETA_C, 'MAX_ELEMS': sp.MAX_ELEMS, 'K_FRAMES': sp.L.K_FRAMES,
       'RC_YU': sp.fdt.RC_YU, 'RC_AH': sp.fdt.RC_AH, 'EPS_AH': sp.fdt.EPS_AH}
for T, I, ph in ((300.0, 0.15, 7.5), (273.0, 0.02, 5.5), (298.0, 0.2, 6.65)):
    out[f'diel{T}'] = sp.fdt.get_dielectric_params(T, I)
params = sp.L.load_residue_params()
out['params'] = params
for ph in (7.5, 5.5):
    q, lam, sig, extra = sp.fdt.calvados_residue_params(SEQ, params, ph=ph)
    out[f'res{ph}'] = (q, lam, sig)
    p = sp.fdt.pair_setup(q, lam, sig)
    out[f'pair{ph}'] = p
    out[f'site{ph}'] = sp.site_matrix(p[0], p[1], q, len(SEQ)).toarray()
for uid, hint in pickle.loads(bytes.fromhex(sys.argv[2])):
    xyz, G, rg = sp.bs.frames_calvados(uid, 300, hint)
    out[f'frames {uid}'] = (xyz, G, rg)
    out[f'terms {uid}'] = sp.calvados_terms(uid, hint)
pickle.dump(out, open(sys.argv[3], 'wb'))
""".replace("SEQ", repr(SEQ))


def same(a, b, path=""):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        return [x for k in a for x in same(a[k], b[k], f"{path}/{k}")]
    if isinstance(a, (tuple, list)):
        assert len(a) == len(b), path
        return [x for i, (u, v) in enumerate(zip(a, b)) for x in same(u, v, f"{path}[{i}]")]
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return [] if np.array_equal(np.asarray(a), np.asarray(b), equal_nan=np.asarray(a).dtype.kind == "f") else [path]
    if hasattr(a, "equals"):
        return [] if a.equals(b) else [path]
    return [] if a == b or (a != a and b != b) else [path]


def main():
    uids = [(u, None) for u in sys.argv[1:]]
    hexa = pickle.dumps(uids).hex()
    res = {}
    for label, path in (("originals", MT), ("helpers", str(ROOT / "helpers"))):
        out = tempfile.mktemp(suffix=".pkl")
        r = subprocess.run([sys.executable, "-B", "-c", CHILD, path, hexa, out], capture_output=True, text=True, cwd=tempfile.gettempdir())
        if r.returncode:
            sys.exit(f"{label} failed:\n{r.stderr[-1500:]}")
        res[label] = pickle.load(open(out, "rb"))
        os.remove(out)
        print(f"{label}: {res[label].pop('file')}")
    diff = same(res["originals"], res["helpers"])
    print(f"compared {len(res['helpers'])} top-level outputs on {len(uids)} trajectories; differences: {len(diff)}")
    print("\n".join(diff[:30]))
    sys.exit(1 if diff else 0)


if __name__ == "__main__":
    main()
