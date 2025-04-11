from typing import Iterable
from random import gauss
import numpy as np
import matplotlib.pyplot as plt

from reactions import *


class Simulation:
    def __init__(self, reactions:Iterable[Reaction]):
        self.reactions = reactions

    def step(self, x:float) -> float:
        ksi = gauss()
        dx_dt = sum(reaction.contribution(x, ksi) for reaction in self.reactions)
        new_x = x + dx_dt
        if new_x < 0: new_x = 0
        return new_x
    
    def run_steps(self, x0:float, n_steps:int) -> Iterable[float]:
        x_t = x0
        xs = [x_t]
        for _ in range(n_steps):
            x_t = self.step(x_t)
            xs.append(x_t)
        return xs
    

if __name__ == "__main__":
    mu = 2
    K = 6
    H = 4.08 
    C = 1
    rho = 1
    dt = 1e-3
    x0 = K
    n_steps = 10000

    _, ax = plt.subplots()
    axs = np.array([[ax]])

    gamma1 = gamma2 = gamma3 = 0
    reactions = [
        LogisticBirth(gamma1, dt, mu),
        LogisticDeath(gamma2, dt, mu, K),
        Harvesting(gamma3, dt, H, C, rho),
    ]
    sim = Simulation(reactions)
    xs = sim.run_steps(x0, n_steps)
    times = [dt * i for i in range(n_steps+1)]
    axs[0,0].plot(times, xs, label=f'gamma1=gamma2=0')

    for gamma in np.logspace(-2, 0, 6):
        gamma1 = gamma2 = gamma
        gamma3 = 0
        reactions = [
            LogisticBirth(gamma1, dt, mu),
            LogisticDeath(gamma2, dt, mu, K),
            Harvesting(gamma3, dt, H, C, rho),
        ]
        sim = Simulation(reactions)
        xs = sim.run_steps(x0, n_steps)
        times = [dt * i for i in range(n_steps+1)]
        axs[0,0].plot(times, xs, label=f'gamma1=gamma2={gamma:.4f}')

    axs[0,0].set_xlabel('Time')
    axs[0,0].set_ylabel('Resources')
    axs[0,0].legend()

    plt.show()