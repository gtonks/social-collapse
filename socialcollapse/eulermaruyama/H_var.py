from typing import Iterable
from math import sqrt
from random import seed, gauss
import numpy as np
import matplotlib.pyplot as plt

from socialcollapse import defaults
from .reactions import *


class Simulation:
    def __init__(self, mu, K, C, rho, dt, H_bar, sigma):
        self.mu_dt = mu * dt
        self.K = K
        self.C_H_bar_dt = C * H_bar * dt
        self.rho = rho
        self.C_sigma_sqrt_dt = C * sigma * sqrt(dt)

    def step(self, x:float) -> float:
        ksi = gauss()
        log_det = logistic_growth(x, self.mu_dt, self.K)
        harvesting_det = harvesting(x, self.C_H_bar_dt, self.rho)
        harvesting_var = harvesting(x, self.C_sigma_sqrt_dt, self.rho) * ksi
        dx = log_det + harvesting_det + harvesting_var
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
    def simulate(mu, K, C, rho, dt, H_bar, sigma, x0, n_steps):
        sim = Simulation(mu, K, C, rho, dt, H_bar, sigma)
        xs = sim.run_steps(x0, n_steps)
        times = [dt * i for i in range(n_steps+1)]
        return xs, times

    

if __name__ == "__main__":
    seed(7)
    mu = defaults.MU
    K = defaults.K
    H_bar = 4.08 
    C = defaults.C
    rho = defaults.RHO
    dt = 1e-3
    x0 = K
    n_steps = 100000

    _, ax = plt.subplots()
    axs = np.array([[ax]])

    sigma = 0
    xs, times = Simulation.simulate(mu, K, C, rho, dt, H_bar, sigma, x0, n_steps)
    axs[0,0].plot(times, xs, label=f'sigma=0')

    for sigma in np.logspace(-2, 0, 6):
        xs, times = Simulation.simulate(mu, K, C, rho, dt, H_bar, sigma, x0, n_steps)
        axs[0,0].plot(times, xs, label=f'{sigma=:.4f}')

    axs[0,0].set_xlabel('Time')
    axs[0,0].set_ylabel('Resources')
    axs[0,0].legend()

    plt.show()