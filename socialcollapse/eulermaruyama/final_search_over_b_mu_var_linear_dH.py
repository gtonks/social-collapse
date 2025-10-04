from time import time
from random import seed
import numpy as np
from .mu_var_linear_dH import Simulation

if __name__ == "__main__":
    seed(24601)
    mu = 1.0
    K = 1.0
    C = 1.031
    rho = 0.04
    r = 2.052
    a = 2.0
    bs = np.concatenate((np.array([11.6, 11.8]), np.arange(12.0, 55.0)))
    dt = 1e-3


    print("Computing deterministic steady states...")

    n_steps = 1_000_000
    sigma = 0
    start = time()

    x_steady = np.empty(bs.size)
    H_steady = np.empty(bs.size)
    for i, b in enumerate(bs):
        xs, Hs, _ = Simulation.simulate(mu, K, C, rho, r, a, b, dt, sigma, 0.9, 0.1, n_steps)
        x_steady[i] = xs[-1]
        H_steady[i] = Hs[-1]
        print(f'\r{i+1}/{bs.size} steady states computed in {time()-start:.2f} seconds', end='')
    print()


    print("Starting trials...")

    n_steps = 10_000
    n_trials = 1000
    sigmas = np.linspace(0, 0.5, 40)

    x_final = np.empty((bs.size, sigmas.size, n_trials))
    H_final = np.empty((bs.size, sigmas.size, n_trials))

    total_trials = n_trials * sigmas.size * bs.size
    current_trial = 0
    start = time()

    for i, b in enumerate(bs):
        for j, sigma in enumerate(sigmas):
            for k in range(n_trials):
                xs, Hs, _ = Simulation.simulate(mu, K, C, rho, r, a, b, dt, sigma, x_steady[i], H_steady[i], n_steps)
                x_final[i, j, k] = xs[-1]
                H_final[i, j, k] = Hs[-1]
                current_trial += 1
            print(f'\r{current_trial}/{total_trials} trials completed in {time()-start:.2f} seconds', end='')
    print()


    file_path = "final_search_over_b.npz"
    print(f"Saving results to {file_path}...")
    np.savez_compressed(file_path,
                        dimensions=np.array(["bs", "sigmas", "trials"]),
                        bs=bs,
                        sigmas=sigmas,
                        x_steady=x_steady,
                        H_steady=H_steady,
                        x_final=x_final,
                        H_final=H_final)