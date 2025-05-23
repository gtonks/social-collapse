from socialcollapse import defaults


def dx_dt(h, x):
    return dx_dt_with_params(h, x, defaults.MU, defaults.K, defaults.C, defaults.RHO)


def dx_dt_with_params(h, x, mu, k, c, rho):
    return mu * x * (1 - x / k) - c * x * h / (rho + x)