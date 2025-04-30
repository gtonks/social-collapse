def logistic_growth(x, mu, K):
    return mu * x * (1 - x / K)


def harvesting(x, HC, rho):
    return -HC * x / (rho + x)