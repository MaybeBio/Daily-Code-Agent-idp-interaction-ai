import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

folder = Path("GA001/curve_1")

def read_numeric(path, min_cols=3):
    rows = []
    for line in path.read_text().splitlines():
        try:
            values = [float(x) for x in line.replace(",", " ").split()]
        except ValueError:
            continue
        if len(values) >= min_cols:
            rows.append(values)

    if not rows:
        raise ValueError(f"No numeric data found in {path}")

    ncols = min(map(len, rows))
    return np.array([row[:ncols] for row in rows])

# 1. SAXS fit
f = next(folder.glob("profiles_*.fit"))
data = read_numeric(f, min_cols=4)

q, exp, err, fit = data[:, :4].T

plt.figure()
plt.errorbar(q, exp, yerr=err, fmt=".", label="Experiment")
plt.plot(q, fit, label="GAJOE fit")
plt.yscale("log")
plt.xlabel(r"$q$ ($\AA^{-1}$)")
plt.ylabel(r"$I(q)$")
plt.legend()
plt.tight_layout()
plt.savefig("saxs_fit.png", dpi=300)

# 2. Structural distributions
for pattern, xlabel, output in [
    ("Rg_distr_*.dat", r"$R_g$ ($\AA$)", "rg.png"),
    ("Size_distr_*.dat", r"$D_{max}$ ($\AA$)", "dmax.png"),
    ("CaCa_distr_*.dat", r"$R_{ee}$ ($\AA$)", "ree.png"),
]:
    files = list(folder.glob(pattern))
    if not files:
        continue

    data = read_numeric(files[0], min_cols=2)

    plt.figure()
    for j in range(1, data.shape[1]):
        plt.plot(data[:, 0], data[:, j], label=f"Column {j+1}")

    plt.xlabel(xlabel)
    plt.ylabel("Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=300)

plt.show()
