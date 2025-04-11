from math import sqrt


class Reaction:
    def __init__(self, gamma, dt):
        self.gamma_sqrt_dt = gamma * sqrt(dt)

    def contribution(self, x:float, ksi:float) -> float:
        raise NotImplementedError
    

class LogisticBirth(Reaction):
    def __init__(self, gamma, dt, mu):
        super().__init__(gamma, dt)
        self.mu = mu

    def contribution(self, x:float, ksi:float) -> float:
        result = self.mu * x * (1 + self.gamma_sqrt_dt * ksi)
        return result
    

class LogisticDeath(Reaction):
    def __init__(self, gamma, dt, mu, K):
        super().__init__(gamma, dt)
        self.mu_over_K = mu / K

    def contribution(self, x:float, ksi:float) -> float:
        result = -self.mu_over_K * x*x * (1 + self.gamma_sqrt_dt * ksi)
        return result
    

class Harvesting(Reaction):
    def __init__(self, gamma, dt, H, C, rho):
        super().__init__(gamma, dt)
        self.HC = H * C
        self.rho = rho

    def contribution(self, x:float, ksi:float) -> float:
        result = -(self.HC * x / (self.rho + x)) * (1 + self.gamma_sqrt_dt * ksi)
        return result