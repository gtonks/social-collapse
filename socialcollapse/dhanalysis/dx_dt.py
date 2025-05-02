from socialcollapse import defaults


def dx_dt(h, x):
    return defaults.MU * x * (1 - x / defaults.K) - (x * h) / (1 + x)