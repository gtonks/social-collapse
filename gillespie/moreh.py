import matplotlib.pyplot as plt

from hs import simulate_plot_hvar


if __name__ == "__main__":
    _, ax = plt.subplots()

    h = 3.1
    simulate_plot_hvar(h, ax, f"H={h}")

    plt.show()