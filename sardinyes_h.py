import sys
import numpy as np
import matplotlib.pyplot as plt

from gillespie3 import run_time
from socialcollapse import get_alpha, get_beta
from reaction import *


if __name__ == "__main__":
    n_samples = 100
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
    for i, H in enumerate(hs):
        print(f'\r{i}/{n_samples}', end='')
        consumption_rxn.alpha = get_alpha(mu, rho, H, C)
        consumption_rxn.beta = get_beta(mu, K, H, C)
        _, A = run_time(A_0, reactions, final_time)
        equilibria[i] = A[-1] / sys_size
    print(f'\r{n_samples}/{n_samples}')
    
    plt.plot(hs, equilibria)
    plt.show()