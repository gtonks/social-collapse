import numpy as np


class Reaction:
    def __init__(self, reactants: int, products: int):
        """
        `reactants`: Number of reactants needed for reaction.
        `products`: Number of products from reaction.
        """
        self.reactants = reactants
        self.products = products

    def propensity(self, A_t: int) -> float:
        """
        Returns the propensity of the reaction at time `t`.

        `A_t`: Number of reactants at time `t`.
        """
        raise NotImplementedError()
    

class LogisticBirth(Reaction):
    def __init__(self):
        super().__init__(0, 1)

    def propensity(self, A_t):
        return 1
    

class LogisticDeath(Reaction):
    def __init__(self, sys_size: float):
        super().__init__(2, 1)
        self.sys_size = sys_size

    def propensity(self, A_t):
        return A_t / self.sys_size
    

class Consumption(Reaction):
    def __init__(self, sys_size: float, alpha: float, beta: float):
        super().__init__(1, 0)
        self.sys_size = sys_size
        self.alpha = alpha
        self.beta = beta

    def propensity(self, A_t):
        if A_t <= 0:
            return 0
        
        denominator = self.alpha + self.beta * A_t
        if denominator == 0 or denominator == np.inf: return 0

        return self.sys_size / denominator    # TODO: Should we multiply by sys_size here, or something else?