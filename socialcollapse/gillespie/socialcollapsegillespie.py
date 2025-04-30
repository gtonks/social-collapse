from math import log
import matplotlib.pyplot as plt
import numpy as np

from reaction import *


class SocialCollapseGillespie:
    def __init__(self, sys_size:float, A:float, B:float):
        """
        Uses the nondimensionalized form dy/dtau = y - y^2 - Ay/(B+y).
        """
        self.reactions = (LogisticBirth(), LogisticDeath(sys_size), Consumption(sys_size, A, B))

    def step(self, y_t:int):
        """
        Computes `tau`, the time to the next reaction, and `A_t_plus_tau`, the number of reactants at time t + tau.

        `y_t`: The number of reactants at time `t`.
        """
        r = np.random.random(size=2)
        propensities = np.array([reaction.propensity(y_t) for reaction in self.reactions])
        alpha0 = propensities.sum()
        if alpha0 <= 0: return np.inf, 0

        tau = -log(r[0]) / alpha0
        y_t_plus_tau = y_t + 1 if r[1] < propensities[0] / alpha0 else y_t - 1
        return tau, y_t_plus_tau

    def run_steps(self, y_0:int, n_steps:int):
        times = [0]
        ys = [y_0]
        for _ in range(n_steps):
            tau, y_t_plus_tau = self.step(ys[-1])
            times.append(times[-1] + tau)
            ys.append(y_t_plus_tau)
        return times, ys

    def run_time(self, y_0:int, final_time:float, return_all:bool=True):
        if return_all:
            times = [0]
            ys = [y_0]
            t, y_t_plus_tau = self.step(y_0)
            while t <= final_time:
                times.append(t)
                ys.append(y_t_plus_tau)
                tau, y_t_plus_tau = self.step(ys[-1])
                t += tau
            times.append(final_time)
            ys.append(ys[-1])
            return times, ys
        
        y_t = y_0
        t, y_t_plus_tau = self.step(y_0)
        while t <= final_time:
            y_t = y_t_plus_tau
            tau, y_t_plus_tau = self.step(y_t)
            t += tau
        return y_t
    

class DimensionalizedSocialCollapseGillespie(SocialCollapseGillespie):
    @staticmethod
    def get_A(H, C, mu, K):
        return H * C / (mu * K)
    
    @staticmethod
    def get_B(rho, K):
        return rho / K
    
    def __init__(self, sys_size, mu, K, C, rho, H):
        super().__init__(sys_size, self.get_A(H, C, mu, K), self.get_B(rho, K))


if __name__ == "__main__":
    sys_size = 500
    y_0 = sys_size
    mu = 2
    K = 6
    H = 5000
    C = 1
    rho = 1
    alg = DimensionalizedSocialCollapseGillespie(sys_size, mu, K, C, rho, H)

    _, (ax1, ax2) = plt.subplots(1, 2)
    for i in range(5):
        np.random.seed(i)
        times, ys = alg.run_steps(y_0, 400000)
        density = np.array(ys) / sys_size
        ax1.plot(times, density, drawstyle='steps-post')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Resource Density')
    ax1.set_xlim(0,2.5)

    for i in range(5):
        np.random.seed(i)
        times, ys = alg.run_time(y_0, 2.5)
        density = np.array(ys) / sys_size
        ax2.plot(times, density, drawstyle='steps-post')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Resource Density')
    plt.tight_layout()
    plt.show()