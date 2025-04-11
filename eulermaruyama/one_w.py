from typing import Iterable
from random import gauss
import matplotlib.pyplot as plt

from reactions import *


class Simulation:
    def __init__(self, reactions:Iterable[Reaction]):
        self.reactions = reactions

    def step(self, x:float) -> float:
        ksi = gauss()
        return sum(reaction.contribution(x, ksi) for reaction in self.reactions)
    
    def run_steps(self, x0:float, n_steps:int) -> Iterable[float]:
        x_t = x0
        xs = [x_t]
        for _ in range(n_steps):
            x_t += self.step(x_t)
            xs.append(x_t)
        return xs
    

if __name__ == "__main__":
    mu = 2
    K = 6
    H = 4.08 
    C = 1
    rho = 1
    dt = 1e-1
    x0 = K
    n_steps = 1000

    gamma1 = 0
    gamma2 = 0
    gamma3 = 0
    reactions = [
        LogisticBirth(gamma1, dt, mu),
        LogisticDeath(gamma2, dt, mu, K),
        Harvesting(gamma3, dt, H, C, rho),
    ]
    sim = Simulation(reactions)
    xs = sim.run_steps(x0, n_steps)
    times = [dt * i for i in range(n_steps+1)]
    plt.plot(times, xs)
    plt.show()