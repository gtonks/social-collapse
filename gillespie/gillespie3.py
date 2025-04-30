from math import log
from collections.abc import Iterable
import matplotlib.pyplot as plt
import numpy as np

from reaction import *


def step(y_t: int, reactions: Iterable[Reaction]):
    """
    Computes `tau`, the time to the next reaction, and `A_t_plus_tau`, the number of reactants at time t + tau.

    `y_t`: The number of reactants at time `t`.
    `reactions`: List of reactions.
    """
    r = np.random.random(size=2)
    propensities = np.array([reaction.propensity(y_t) for reaction in reactions])
    # print(f"{A_t=}\n{propensities=}")
    alpha0 = propensities.sum()
    if alpha0 <= 0: return np.inf, 0

    tau = -log(r[0]) / alpha0

    propensity_ratios = propensities / alpha0
    j = 1
    sum_alpha1_alphaj = propensity_ratios[0]
    while sum_alpha1_alphaj < r[1]:
        j += 1
        sum_alpha1_alphaj += propensity_ratios[j-1]
    reaction = reactions[j-1]
    A_t_plus_tau = y_t - reaction.reactants + reaction.products
    return tau, A_t_plus_tau


def run_steps(A_0: int, reactions: Iterable[Reaction], n_steps: int):
    times = [0]
    A = [A_0]
    for _ in range(n_steps):
        tau, A_t_plus_tau = step(A[-1], reactions)
        times.append(times[-1] + tau)
        A.append(A_t_plus_tau)
    return times, A


def run_time(A_0: int, reactions: Iterable[Reaction], final_time: float, return_all:bool=True):
    if return_all:
        times = [0]
        A = [A_0]
        current_time, A_t_plus_tau = step(A_0, reactions)
        while current_time <= final_time:
            times.append(current_time)
            A.append(A_t_plus_tau)
            tau, A_t_plus_tau = step(A[-1], reactions)
            current_time += tau
        return times, A
    
    time = 0
    A_t = A_0
    current_time, A_t_plus_tau = step(A_0, reactions)
    while current_time <= final_time:
        time = current_time
        A_t = A_t_plus_tau
        tau, A_t_plus_tau = step(A_t, reactions)
        current_time += tau
    return A_t


if __name__ == "__main__":
    sys_size = 500
    A_0 = sys_size
    mu = 2
    K = 6
    H = 5
    C = 1
    rho = 1
    alpha = mu * rho / (H * C)
    beta = mu * K / (H * C)
    print((alpha + beta)**2 >= 4*beta)
    reactions = [
        LogisticBirth(),
        LogisticDeath(sys_size),
        Consumption(sys_size, alpha, beta)
    ]

    _, (ax1, ax2) = plt.subplots(1, 2)
    for i in range(5):
        np.random.seed(i)
        times, A = run_steps(A_0, reactions, 40000)
        density = np.array(A) / sys_size
        ax1.plot(times, density, drawstyle='steps-post')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Resource Density')
    # for i in range(5):
    #     np.random.seed(i)
    #     times, A = run_time(A_0, reactions, 20000)
    #     density = np.array(A) / sys_size
    #     ax2.plot(times, density, drawstyle='steps-post')
    # ax2.set_xlabel('Time')
    # ax2.set_ylabel('Resource Density')
    plt.show()