import numpy as np
import matplotlib.pyplot as plt

from .dx_dt import dx_dt
from . import vectors


def dh_dt(h, x):
    return np.zeros_like(x)

if __name__ == "__main__":
    _, axs = plt.subplots(1,2)
    max_h = 5
    max_x = 6
    
    ax = axs[0]
    vectors.plot_field(dh_dt, dx_dt, ax, max_h, max_x)
    ax.set_xlabel('h')
    ax.set_ylabel('x')

    ax = axs[1]
    vectors.plot_streamlines(dh_dt, dx_dt, ax, max_h, max_x)
    ax.set_xlabel('h')
    ax.set_ylabel('x')

    plt.show()