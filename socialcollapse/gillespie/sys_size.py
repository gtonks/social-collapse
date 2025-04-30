import matplotlib.pyplot as plt

from socialcollapse import simulate_plot
from reaction import *


def simulate_plot_syssize(sys_size, time, ax, title=""):
    simulate_plot(ax, time, sys_size, A_0=sys_size, mu=2, K=6, H=3.2, C=1, rho=1, title=title)
    

if __name__ == "__main__":
    _, axs = plt.subplots(1, 3)

    simulate_plot_syssize(sys_size=5e2, time=1e4, ax=axs[0], title="System size = 500")
    simulate_plot_syssize(sys_size=5e3, time=1e5, ax=axs[1], title="System size = 5000")
    simulate_plot_syssize(sys_size=5e4, time=1e6, ax=axs[2], title="System size = 50000")

    plt.show()