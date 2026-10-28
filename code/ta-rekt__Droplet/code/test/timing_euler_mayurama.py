import numpy as np

def euler_maruyama(x0, kappa, gamma, D, dt, n, seed=None):
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = x0
    xi = rng.standard_normal(n - 1)
    for i in range(n - 1):
        x[i+1] = x[i] - (kappa/gamma)*x[i]*dt + np.sqrt(2*D*dt)*xi[i]
    return x

if __name__ == "__main__":

    import time

    start = time.perf_counter()
    
    x0 = 0.0
    kappa = 1.0
    gamma = 1.0
    D = 1.0
    dt = 0.01 * (gamma / kappa)
    n = 10_000_000

    euler_maruyama(x0, kappa, gamma, D, dt, n, seed=None)

    end = time.perf_counter()

    print(end - start)