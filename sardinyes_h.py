import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib.pyplot as plt

from gillespie3 import run_time
from socialcollapse import get_alpha, get_beta
from reaction import *


if __name__ == "__main__":
    _, ax = plt.subplots()

    n_samples = 10
    hs = np.linspace(0, 5, n_samples)
    equilibria = np.zeros(n_samples)
    sys_size = 500
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

    start = time.time()
    for i, H in enumerate(hs):
        print(f'\r{i}/{n_samples}', end='')
        consumption_rxn.alpha = get_alpha(mu, rho, H, C)
        consumption_rxn.beta = get_beta(mu, K, H, C)
        A_final = run_time(A_0, reactions, final_time, return_all=False)
        equilibria[i] = A_final / sys_size
    print(f'\r{n_samples}/{n_samples}')
    end=time.time()
    print(f'Total time: {end-start}')
    
    ax.plot(hs, equilibria)
    ax.set_xlim(hs.min(), hs.max())
    ax.set_xlabel('H')
    ax.set_ylabel('Resource Density at Equilibrium')

    plt.show()