import numpy as np

from gillespie3 import run_time
from reaction import *


def get_alpha(mu, rho, H, C):
    denominator = H * C
    return mu * rho / denominator if denominator != 0 else np.inf


def get_beta(mu, K, H, C):
    denominator = H * C
    return mu * K / denominator if denominator != 0 else np.inf


def simulate_plot(ax, time=1e6, sys_size=5e4, A_0=5e4, mu=2, K=6, H=4.08, C=1, rho=1, title=""):
    alpha = get_alpha(mu, rho, H, C)
    beta = get_beta(mu, K, H, C)
    print((alpha + beta)**2 >= 4*beta)
    reactions = [
        LogisticBirth(),
        LogisticDeath(sys_size),
        Consumption(sys_size, alpha, beta)
    ]

    for i in range(5):
        np.random.seed(i)
        times, A = run_time(A_0, reactions, time)
        print(f"{len(A)=}")
        density = np.array(A) / sys_size
        ax.plot(times, density, drawstyle='steps-post')
    ax.set_title(title)
    ax.set_xlim(0, time)
    ax.set_xlabel('Time')
    ax.set_ylim(0, 1)
    ax.set_ylabel('Resource Density')