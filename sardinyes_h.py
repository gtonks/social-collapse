import time
import hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib.pyplot as plt

from gillespie3 import run_time
from socialcollapse import get_alpha, get_beta
from reaction import *


def seed_H(H):
    np.random.seed(int(hashlib.md5(str(H).encode()).hexdigest(), 16) % 2**32)


def run_series(hs):
    n_samples = len(hs)
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
        seed_H(H)
        print(f'\r{i}/{n_samples}', end='')
        consumption_rxn.alpha = get_alpha(mu, rho, H, C)
        consumption_rxn.beta = get_beta(mu, K, H, C)
        A_final = run_time(A_0, reactions, final_time, return_all=False)
        equilibria[i] = A_final / sys_size
    print(f'\r{n_samples}/{n_samples}')
    return equilibria


def get_final(H):
    seed_H(H)
    sys_size = 500
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


def run_parallel(hs):
    n_samples = len(hs)
    equilibria = np.zeros(n_samples)
    with ProcessPoolExecutor() as executor:
        equilibria = list(executor.map(get_final, hs))
    return equilibria


if __name__ == "__main__":
    pass