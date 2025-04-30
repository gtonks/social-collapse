from typing import Iterable
from math import sqrt
from random import seed, gauss
import numpy as np
import matplotlib.pyplot as plt

from socialcollapse import defaults
from .reactions import *


def logistic_growth(x, mu, K):
    return mu * x * (1 - x / K)


def harvesting(x, HC, rho):
    return -HC * x / (rho + x)


class Simulation:
    def __init__(self, K, H, C, rho, dt, mu_bar, sigma):
        self.mu_bar_dt = mu_bar * dt
        self.K = K
        self.H_C_dt = H * C * dt
        self.rho = rho
        self.sigma_sqrt_dt = sigma * sqrt(dt)

    def step(self, x:float) -> float:
        ksi = gauss()
        log_det = logistic_growth(x, self.mu_bar_dt, self.K)
        harvesting_det = harvesting(x, self.H_C_dt, self.rho)
        log_var = logistic_growth(x, self.sigma_sqrt_dt, self.K) * ksi
        dx = log_det + harvesting_det + log_var
        new_x = x + dx
        if new_x < 0: new_x = 0
        return new_x

    def run_steps(self, x0:float, n_steps:int) -> Iterable[float]:
        x_t = x0
        xs = [x_t]
        for _ in range(n_steps):
            x_t = self.step(x_t)
            xs.append(x_t)
        return xs

    @staticmethod
    def simulate(K, H, C, rho, dt, mu_bar, sigma, x0, n_steps):
        sim = Simulation(K, H, C, rho, dt, mu_bar, sigma)
        xs = sim.run_steps(x0, n_steps)
        times = [dt * i for i in range(n_steps+1)]
        return xs, times

    

if __name__ == "__main__":
    seed(7)
    mu_bar = defaults.MU
    K = defaults.K
    H = 4.08 
    C = defaults.C
    rho = defaults.RHO
    dt = 1e-3
    x0 = K
    n_steps = 100000

    _, ax = plt.subplots()
    axs = np.array([[ax]])

    sigma = 0.01
    xs, times = Simulation.simulate(K, H, C, rho, dt, mu_bar, sigma, x0, n_steps)
    axs[0,0].plot(times, xs, label=f'sigma=0')

    for sigma in np.logspace(-2, 0, 6):
        xs, times = Simulation.simulate(K, H, C, rho, dt, mu_bar, sigma, x0, n_steps)
        axs[0,0].plot(times, xs, label=f'{sigma=:.4f}')

    axs[0,0].set_xlabel('Time')
    axs[0,0].set_ylabel('Resources')
    axs[0,0].legend()

    plt.show()