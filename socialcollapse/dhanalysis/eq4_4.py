from math import log10
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec

from socialcollapse import defaults
from .dx_dt import dx_dt_with_params
from . import vectors


def gamma_with_params(x, alpha, beta, k):
    return beta - (beta - alpha) * x / k


def dh_dt_with_params(h, x, r, alpha, beta, k):
    return h * (r - gamma_with_params(x, alpha, beta, k) * h)


if __name__ == "__main__":
    init_r = 1
    init_alpha = 0.001
    init_beta = 1
    init_mu = defaults.MU
    init_k = defaults.K
    init_c = defaults.C
    init_rho = defaults.RHO

    max_h = 10
    max_x = 8

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    # Initial vector field
    dh_dt = lambda h, x: dh_dt_with_params(h, x, init_r, init_alpha, init_beta, init_k)
    dx_dt = lambda h, x: dx_dt_with_params(h, x, init_mu, init_k, init_c, init_rho)

    stream, sx_grid, sy_grid = vectors.plot_streamlines(dh_dt, dx_dt, ax, max_h, max_x)
    ax.set_xlabel('h')
    ax.set_ylabel('x')
    ax.set_title(f'Stream: mu={init_mu:.2f}, r={init_r:.2f}, alpha={init_alpha:.2f}, beta={init_beta:.2f}, c={init_c:.2f}, k={init_k:.2f}')
    
    # Slider axes
    ax_grid = gridspec.GridSpec(2, 3, left=0.1, bottom=0.05, right=0.9, top=0.15, wspace=0.4, hspace=0.4)
    ax_mu_slider = plt.subplot(ax_grid[0,0])
    ax_r_slider = plt.subplot(ax_grid[1,0])
    ax_alpha_slider = plt.subplot(ax_grid[0,1])
    ax_beta_slider = plt.subplot(ax_grid[1,1])
    ax_c_slider = plt.subplot(ax_grid[0,2])
    ax_k_slider = plt.subplot(ax_grid[1,2])

    mu_slider = Slider(ax_mu_slider, 'mu (log)', -2, 2, valinit=log10(init_mu))
    r_slider = Slider(ax_r_slider, 'r (log)', -2, 2, valinit=log10(init_r))
    alpha_slider = Slider(ax_alpha_slider, 'alpha (log)', -3, 0, valinit=log10(init_alpha))
    beta_slider = Slider(ax_beta_slider, 'beta (log)', 0, 3, valinit=log10(init_beta))
    c_slider = Slider(ax_c_slider, 'c (log)', -1, 1, valinit=log10(init_c))
    k_slider = Slider(ax_k_slider, 'k (log)', 0, 1, valinit=log10(init_k))

    def update(val):
        mu = 10 ** mu_slider.val
        r = 10 ** r_slider.val
        alpha = 10 ** alpha_slider.val
        beta = 10 ** beta_slider.val
        c = 10 ** c_slider.val
        k = 10 ** k_slider.val
        rho = init_rho

        ax.cla()
        new_dh_s = dh_dt_with_params(sx_grid, sy_grid, r, alpha, beta, k)
        new_dx_s = dx_dt_with_params(sx_grid, sy_grid, mu, k, c, rho)
        ax.streamplot(sx_grid, sy_grid, new_dh_s, new_dx_s)
        ax.set_xlabel('h')
        ax.set_ylabel('x')
        ax.set_title(f'Stream: {mu=:.2f}, {r=:.2f}, {alpha=:.2f}, {beta=:.2f}, {c=:.2f}, {k=:.2f}')

        fig.canvas.draw_idle()

    mu_slider.on_changed(update)
    r_slider.on_changed(update)
    alpha_slider.on_changed(update)
    beta_slider.on_changed(update)
    c_slider.on_changed(update)
    k_slider.on_changed(update)

    plt.show()
