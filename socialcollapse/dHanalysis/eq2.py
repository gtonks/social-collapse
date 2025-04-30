import numpy as np
import matplotlib.pyplot as plt

from .dx_dt import dx_dt
from . import vectors


def dh_dt_with_params(h, x, mu_h, k_h):
    return mu_h * x * h * (1 - h * h / k_h)

if __name__ == "__main__":
    mu_hs = [0.1, 1, 10]
    k_h = 4

    _, axs = plt.subplots(len(mu_hs),2)
    max_h = 5
    max_x = 6
    
    for i, mu_h in enumerate(mu_hs):
        dh_dt = lambda h, x: dh_dt_with_params(h, x, mu_h, k_h)

        ax = axs[i, 0]
        vectors.plot_field(dh_dt, dx_dt, ax, max_h, max_x)
        ax.set_xlabel('h')
        ax.set_ylabel('x')
        ax.set_title(f'{mu_h=}')

        ax = axs[i, 1]
        vectors.plot_streamlines(dh_dt, dx_dt, ax, max_h, max_x)
        ax.set_xlabel('h')
        ax.set_ylabel('x')
        ax.set_title(f'{mu_h=}')

    plt.tight_layout()
    plt.show()