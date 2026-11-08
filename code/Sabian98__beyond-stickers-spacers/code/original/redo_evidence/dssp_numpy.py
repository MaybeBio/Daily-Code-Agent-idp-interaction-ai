"""
dssp_numpy.py — secondary structure by the DSSP algorithm of Kabsch & Sander (Biopolymers 1983), in numpy, from the
backbone atoms N, CA, C, O of every model of a multi-model PDB file (no DSSP program on this Mac).

  Hydrogen bond  electrostatic energy E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol < -0.5, with the
                 amide H placed as in DSSP (N + unit vector of C(i-1) -> O(i-1)); no H on proline or after a chain break
                 (C(i-1)-N(i) > 2.5 A). Deposited hydrogens are ignored, as in DSSP.
  n-turn(i)      H-bond from CO(i) to NH(i+n), n = 3, 4, 5.
  Helices        two consecutive n-turns at i-1 and i: residues i..i+n-1 get H (n = 4), G (n = 3) or I (n = 5).
  Bridges        parallel:      [Hb(i-1, j) and Hb(j, i+1)] or [Hb(j-1, i) and Hb(i, j+1)]
                 antiparallel:  [Hb(i, j) and Hb(j, i)] or [Hb(i-1, j+1) and Hb(j-1, i+1)],   |i - j| > 2,
                 where Hb(a, b) = H-bond from CO(a) to NH(b). A bridge with a consecutive partner bridge of the same
                 type (i+-1 with j+-1) is a ladder residue E, otherwise an isolated bridge B (bulges not merged).
  Priority       H > B > E > G > I (T and S not assigned).
Self-test: python3 dssp_numpy.py  (ideal alpha helix, ideal strand)
"""
import numpy as np

Q = 0.084 * 332.0


def read_backbone(path, n_expected=None, max_models=None):
    """List of models, each an array (n_res, 4, 3) of N, CA, C, O (NaN where missing), plus the residue names."""
    models, cur, order, names = [], None, [], []

    def close():
        if cur:
            keys = [k for k in order if 'CA' in cur[k]]
            arr = np.full((len(keys), 4, 3), np.nan)
            for r, k in enumerate(keys):
                for a, nm in enumerate(('N', 'CA', 'C', 'O')):
                    if nm in cur[k]:
                        arr[r, a] = cur[k][nm]
                    elif nm == 'O' and 'OXT' in cur[k]:
                        arr[r, a] = cur[k]['OXT']
            models.append(arr)
            if not names:
                names.extend(cur[k]['_rn'] for k in keys)

    for line in open(path):
        if line.startswith('MODEL'):
            close()
            cur, order = {}, []
        elif line.startswith(('ATOM', 'HETATM')) and line[16] in ' A':
            if cur is None:
                cur, order = {}, []
            k = (line[21], line[22:27])
            if k not in cur:
                cur[k] = {'_rn': line[17:20]}
                order.append(k)
            cur[k][line[12:16].strip()] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        elif line.startswith('ENDMDL'):
            close()
            cur, order = None, []
    close()
    if n_expected is not None:
        models = [m for m in models if len(m) == n_expected]
    if max_models and len(models) > max_models:
        models = models[::int(np.ceil(len(models) / max_models))]
    return models, names


def hbond_matrix(bb, names):
    """hb[a, b] = True if CO of residue a is H-bonded to NH of residue b."""
    N, C, O = bb[:, 0], bb[:, 2], bb[:, 3]
    n = len(bb)
    Hpos = np.full((n, 3), np.nan)
    v = C[:-1] - O[:-1]
    v = v / np.linalg.norm(v, axis=1, keepdims=True)
    Hpos[1:] = N[1:] + v
    brk = np.linalg.norm(N[1:] - C[:-1], axis=1) > 2.5
    Hpos[1:][brk] = np.nan
    pro = np.array([r == 'PRO' for r in names])
    Hpos[pro] = np.nan

    def d(x, y):
        return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        E = Q * (1 / d(O, N) + 1 / d(C, Hpos) - 1 / d(O, Hpos) - 1 / d(C, N))   # [acceptor a, donor b]
    hb = E < -0.5
    hb[np.isnan(E)] = False
    idx = np.arange(n)
    hb[np.abs(idx[:, None] - idx[None, :]) < 2] = False
    return hb


def assign(bb, names):
    """One-letter DSSP codes (H, B, E, G, I or '-') for one model."""
    n = len(bb)
    hb = hbond_matrix(bb, names)
    ss = np.array(['-'] * n)
    helix = {}
    for nn, code in ((4, 'H'), (3, 'G'), (5, 'I')):
        turn = np.zeros(n, bool)
        turn[:n - nn] = hb[np.arange(n - nn), np.arange(nn, n)]
        h = np.zeros(n, bool)
        for i in range(1, n - nn + 1):
            if turn[i - 1] and turn[i]:
                h[i:i + nn] = True
        helix[code] = h

    def Hb(a, b):
        return 0 <= a < n and 0 <= b < n and hb[a, b]
    par, anti = set(), set()
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if abs(i - j) <= 2:
                continue
            if (Hb(i - 1, j) and Hb(j, i + 1)) or (Hb(j - 1, i) and Hb(i, j + 1)):
                par.add((i, j))
            if (Hb(i, j) and Hb(j, i)) or (Hb(i - 1, j + 1) and Hb(j - 1, i + 1)):
                anti.add((i, j))
    bridge, ladder = np.zeros(n, bool), np.zeros(n, bool)
    for (i, j) in par:
        bridge[i] = True
        if (i + 1, j + 1) in par or (i - 1, j - 1) in par:
            ladder[i] = True
    for (i, j) in anti:
        bridge[i] = True
        if (i + 1, j - 1) in anti or (i - 1, j + 1) in anti:
            ladder[i] = True
    for code, mask in (('I', helix['I']), ('G', helix['G']), ('E', ladder), ('B', bridge & ~ladder), ('H', helix['H'])):
        ss[mask] = code
    return ''.join(ss)


def ensemble_fractions(path, n_expected=None, max_models=1000):
    """Per-residue fraction of models in helix (H, G, I), alpha helix (H) and strand (E, B); number of models."""
    models, names = read_backbone(path, n_expected, max_models)
    if not models:
        return None
    S = np.array([list(assign(m, names)) for m in models])
    return (np.isin(S, list('HGI')).mean(0), (S == 'H').mean(0), np.isin(S, list('EB')).mean(0), len(models), S)


def _build(phi, psi, n):
    """Ideal backbone by NeRF (standard bond lengths and angles), residue names ALA."""
    def place(a, b, c, bond, ang, tor):
        ang, tor = np.radians(ang), np.radians(tor)
        bc = c - b
        bc /= np.linalg.norm(bc)
        nv = np.cross(b - a, bc)
        nv /= np.linalg.norm(nv)
        m = np.array([bc, np.cross(nv, bc), nv]).T
        d2 = np.array([-bond * np.cos(ang), bond * np.sin(ang) * np.cos(tor), bond * np.sin(ang) * np.sin(tor)])
        return c + m @ d2
    N, CA, C = np.array([0., 0, 0]), np.array([1.458, 0, 0]), None
    C = CA + 1.525 * np.array([np.cos(np.radians(180 - 111.2)), np.sin(np.radians(180 - 111.2)), 0])
    atoms = [[N, CA, C]]
    for i in range(1, n):
        pN, pCA, pC = atoms[-1]
        Nn = place(pN, pCA, pC, 1.329, 116.2, psi)
        CAn = place(pCA, pC, Nn, 1.458, 121.7, 180.0)
        Cn = place(pC, Nn, CAn, 1.525, 111.2, phi)
        atoms.append([Nn, CAn, Cn])
    bb = np.full((n, 4, 3), np.nan)
    for i in range(n):
        Nn, CAn, Cn = atoms[i]
        nxtN = atoms[i + 1][0] if i + 1 < n else place(Nn, CAn, Cn, 1.329, 116.2, psi)
        bb[i, :3] = atoms[i]
        bb[i, 3] = place(nxtN, CAn, Cn, 1.231, 120.5, 180.0)     # O trans to N(i+1) about CA-C
    return bb, ['ALA'] * n


if __name__ == '__main__':
    bb, nm = _build(-57.0, -47.0, 20)
    print('ideal alpha helix  ', assign(bb, nm))
    bb, nm = _build(-120.0, 130.0, 20)
    print('ideal single strand', assign(bb, nm))
