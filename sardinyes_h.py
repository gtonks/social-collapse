import time
import hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib.pyplot as plt

from gillespie3 import run_time
from socialcollapse import get_alpha, get_beta
from reaction import *


def seed(H, sys_size):
    np.random.seed(int(hashlib.md5(str(H+sys_size).encode()).hexdigest(), 16) % 2**32)


def run_series(hs, sys_size):
    n_samples = len(hs)
    equilibria = np.zeros(n_samples)
    final_time = 1e6
    A_0 = sys_size
    mu = 2
    K = 6
    C = 1
    rho = 1
    consumption_rxn = Consumption(sys_size, 0, 0)
    reactions = [
        LogisticBirth(),
        LogisticDeath(sys_size),
        consumption_rxn
    ]

    for i, H in enumerate(hs):
        seed(H, sys_size)
        print(f'\r{i}/{n_samples}', end='')
        consumption_rxn.alpha = get_alpha(mu, rho, H, C)
        consumption_rxn.beta = get_beta(mu, K, H, C)
        A_final = run_time(A_0, reactions, final_time, return_all=False)
        equilibria[i] = A_final / sys_size
    print(f'\r{n_samples}/{n_samples}')
    return equilibria


def get_final(H, sys_size):
    seed(H, sys_size)
    final_time = 1e6
    A_0 = sys_size
    mu = 2
    K = 6
    C = 1
    rho = 1
    alpha = get_alpha(mu, rho, H, C)
    beta = get_beta(mu, K, H, C)
    reactions = [
        LogisticBirth(),
        LogisticDeath(sys_size),
        Consumption(sys_size, alpha, beta)
    ]
    A_final = run_time(A_0, reactions, final_time, return_all=False)
    return A_final / sys_size


def run_parallel(hs, sys_size, print_updates=None):
    n_samples = len(hs)
    equilibria = [0] * n_samples

    if print_updates:
        h_batches = np.array_split(hs, print_updates)
        with ProcessPoolExecutor() as executor:
            n_completed = 0
            for i, batch in enumerate(h_batches):
                print(f'\r{100 * i / print_updates:.1f}%', end='')
                batch_size = len(batch)
                equilibria[n_completed:n_completed+batch_size] = executor.map(get_final, batch, [sys_size]*batch_size)
                n_completed += batch_size
        print(f'\r{100:.1f}%')

    else:
        with ProcessPoolExecutor() as executor:
            equilibria = list(executor.map(get_final, hs, [sys_size]*n_samples))
    return equilibria


if __name__ == "__main__":
    _, ax = plt.subplots()

    n_samples = 100
    hs = np.linspace(0, 5, n_samples)
    sys_sizes = [500, 5_000, 50_000]

    for sys_size in sys_sizes:
        print(f'Begin {sys_size=}')
        start = time.time()
        equilibria = run_parallel(hs, sys_size, print_updates=10)
        end=time.time()
        print(f'Simulation time: {end-start}')        
        ax.plot(hs, equilibria, 'o', label=f'{sys_size=}')

    ax.set_xlim(hs.min(), hs.max())
    ax.set_xlabel('H')
    ax.set_ylabel('Resource Density at Equilibrium')
    ax.legend()

    plt.show()