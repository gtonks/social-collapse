import matplotlib.pyplot as plt

from socialcollapse import simulate_plot
from reaction import *


def simulate_plot_hvar(h, ax, title=""):
    simulate_plot(ax, time=1e6, sys_size=5e4, A_0=5e4, mu=2, K=6, H=h, C=1, rho=1, title=title)


if __name__ == "__main__":
    _, axs = plt.subplots(2, 2)

    for i, h in enumerate([2.8, 3.0, 3.2, 3.4]):
        simulate_plot_hvar(h, axs[i//2, i%2], f"H={h}")

    plt.tight_layout()
    plt.show()