import numpy as np


class Reaction:
    def __init__(self, reactants: int, products: int):
        """
        `reactants`: Number of reactants needed for reaction.
        `products`: Number of products from reaction.
        """
        self.reactants = reactants
        self.products = products

    def propensity(self, y_t: int) -> float:
        """
        Returns the propensity of the reaction at time `t`.

        `y_t`: Number of reactants at time `t`.
        """
        raise NotImplementedError()
    

class LogisticBirth(Reaction):
    def __init__(self):
        super().__init__(0, 1)

    def propensity(self, y_t):
        return y_t
    

class LogisticDeath(Reaction):
    def __init__(self, sys_size: float):
        super().__init__(2, 1)
        self.two_sys_size = 2 * sys_size

    def propensity(self, y_t):
        return y_t * (y_t - 1) / self.two_sys_size
    

class Consumption(Reaction):
    def __init__(self, sys_size: float, A: float, B: float):
        super().__init__(1, 0)
        self.sys_size = sys_size
        self.two_sys_size_inv = 1 / (2 * sys_size)
        self.A = A
        self.B = B

    def propensity(self, y_t):
        if y_t <= 0:
            return 0
        
        denominator = self.B + y_t
        if denominator == 0 or denominator == np.inf: return 0

        # return y_t - y_t*(y_t-1)*self.two_sys_size_inv
        # return y_t * (y_t - y_t * self.two_sys_size_inv + self.two_sys_size_inv + self.A / denominator)
        return self.A * y_t / denominator    # TODO: Should we multiply by sys_size here, or something else?