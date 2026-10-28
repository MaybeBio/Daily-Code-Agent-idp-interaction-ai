"""Minimal reproduction of the Rouse MSD fitting failure."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from scipy import special
from viz_util.fit_model import fit

TRUE = {"gamma": 0.5, "J": 30.0}

def rouse_msd2(t, params):
    gamma, J = params["gamma"], params["J"]
    tau = (J / gamma) ** 2
    return (2*gamma*t**(1/2) * (1 - np.exp(-tau/(np.pi*t)))
            + 2*J * special.erfc((tau/(np.pi*t))**(1/2)))

time = np.linspace(1, 1000, 200)
x = time[:, np.newaxis]             # (n, 1)  — as used in notebook
y_true = rouse_msd2(x, TRUE)

print("x shape   :", x.shape)
print("y_true shape:", y_true.shape)

# What does model(x[0], p) return?
sample_out = rouse_msd2(x[0], TRUE)
print("model(x[0]) shape:", np.asarray(sample_out).shape, "value:", sample_out)

# What does the internal _f produce?
n = len(x)
p_dict = TRUE
loop_out = np.array([rouse_msd2(x[i], p_dict) for i in range(n)])
print("loop output shape:", loop_out.shape)   # expect (n,) but likely (n,1)
