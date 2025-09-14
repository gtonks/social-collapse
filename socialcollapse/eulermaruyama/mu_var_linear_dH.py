from typing import Iterable
from math import sqrt
from random import seed, gauss
import numpy as np
import matplotlib.pyplot as plt

from .reactions import *


def gamma(x, b, b_minus_a_over_K):
    return b - b_minus_a_over_K * x


def dH_dt(x, H, r, b, b_minus_a_over_K):
    return (r - gamma(x, b, b_minus_a_over_K) * H) * H


class Simulation:
    def __init__(self, mu, K, C, rho, r, a, b, dt, sigma):
        self.dt = dt
        self.mu_dt = mu * dt
        self.K = K
        self.C_dt = C * dt
        self.rho = rho
        self.C_sigma_sqrt_dt = C * sigma * sqrt(dt)
        self.r = r
        self.b = b
        self.b_minus_a_over_K = (b - a) / K

    def step_x(self, x:float, H:float) -> float:
        ksi = gauss(0, 1)
        log_det = logistic_growth(x, self.mu_dt, self.K)
        harvesting_det = harvesting(x, self.C_dt * H, self.rho)
        harvesting_var = harvesting(x, self.C_sigma_sqrt_dt, self.rho) * ksi
        dx = log_det + harvesting_det + harvesting_var
        new_x = x + dx
        if new_x < 0: new_x = 0
        return new_x
    
    def step_H(self, x:float, H:float) -> float:
        dH = dH_dt(x, H, self.r, self.b, self.b_minus_a_over_K) * self.dt
        new_H = H + dH
        if new_H < 0: new_H = 0
        return new_H

    def run_steps(self, x0:float, H0:float, n_steps:int) -> Iterable[float]:
        x_t = x0
        H_t = H0
        xs = [x_t]
        Hs = [H_t]
        for _ in range(n_steps):
            x_t = self.step_x(x_t, H_t)
            H_t = self.step_H(x_t, H_t)
            xs.append(x_t)
            Hs.append(H_t)
        return xs, Hs

    @staticmethod
    def simulate(mu, K, C, rho, r, a, b, dt, sigma, x0, H0, n_steps):
        sim = Simulation(mu, K, C, rho, r, a, b, dt, sigma)
        xs, Hs = sim.run_steps(x0, H0, n_steps)
        times = [dt * i for i in range(n_steps+1)]
        return xs, Hs, times

    

if __name__ == "__main__":
    seed(7)
    mu = 1.0
    K = 1.0
    C = 1.031
    rho = 0.04
    r = 2.052
    a = 2.0
    b = 11.8
    dt = 1e-3
    x0 = 0.9
    H0 = 0.1
    n_steps = 10000

    _, axs = plt.subplots(1,3)
    axs = axs.reshape((1,3))

    sigma = 0
    xs, Hs, times = Simulation.simulate(mu, K, C, rho, r, a, b, dt, sigma, x0, H0, n_steps)
    print(f'Final deterministic state: x={xs[-1]}, H={Hs[-1]}')
    axs[0,0].plot(times, xs, label=f'sigma=0')
    axs[0,1].plot(times, Hs, label=f'sigma=0')
    axs[0,2].plot(xs, Hs, label=f'sigma=0')

    for sigma in np.logspace(-2, 0, 6):
        xs, Hs, times = Simulation.simulate(mu, K, C, rho, r, a, b, dt, sigma, x0, H0, n_steps)
        axs[0,0].plot(times, xs, label=f'{sigma=:.4f}')
        axs[0,1].plot(times, Hs, label=f'{sigma=:.4f}')
        axs[0,2].plot(xs, Hs, label=f'{sigma=:.4f}')

    axs[0,0].set_xlabel('Time')
    axs[0,0].set_ylabel('Resources')
    axs[0,0].legend()

    axs[0,1].set_xlabel('Time')
    axs[0,1].set_ylabel('Humans')
    axs[0,1].legend()

    axs[0,2].set_xlabel('Resources')
    axs[0,2].set_ylabel('Humans')
    axs[0,2].legend()

    plt.show()