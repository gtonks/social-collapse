import numpy as np
import matplotlib.pyplot as plt

from socialcollapse import defaults
from . import vectors

def dx_dt(h, x):
    return defaults.MU * x * (1 - x / defaults.K) - (x * h) / (1 + x)

def dh_dt(h, x):
    return np.zeros_like(x)

if __name__ == "__main__":
    _, axs = plt.subplots(1,2)
    
    ax = axs[0]
    vectors.plot_field(dh_dt, dx_dt, ax)
    ax.set_xlabel('h')
    ax.set_ylabel('x')

    ax = axs[1]
    vectors.plot_streamlines(dh_dt, dx_dt, ax)
    ax.set_xlabel('h')
    ax.set_ylabel('x')

    plt.show()