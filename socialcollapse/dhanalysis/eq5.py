from math import log10
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec

from socialcollapse import defaults
from .dx_dt import dx_dt_with_params
from . import vectors


def gamma_with_params(x, alpha, beta, k):
    return beta - (beta - alpha) * x / k


def dh_dt_with_params(h, x, b, rho, d, alpha, beta, k):
    return h * (b * x / (rho + x) - d - gamma_with_params(x, alpha, beta, k) * h)


def find_fixed_points(params, guesses):
    b, d, alpha, beta, k, mu, c, rho = params
    def system(vars):
        h, x = vars
        return [
            dh_dt_with_params(h, x, b, rho, d, alpha, beta, k),
            dx_dt_with_params(h, x, mu, k, c, rho)
        ]
    results = []
    for guess in guesses:
        sol = root(system, guess)
        if sol.success:
            h_star, x_star = sol.x
            results.append((h_star, x_star))
    return results


def jacobian(h, x, params):
    b, d, alpha, beta, k, mu, c, rho = params
    # Partial derivatives
    gamma = gamma_with_params(x, alpha, beta, k)
    dgamma_dx = -(beta - alpha) / k

    # dh/dt partials
    # ∂(dh/dt)/∂h
    d_dh = b * x / (rho + x) - d - 2 * gamma * h
    # ∂(dh/dt)/∂x
    d_dx = h * (b * rho / (rho + x) ** 2 - dgamma_dx * h)

    # dx/dt partials
    # ∂(dx/dt)/∂h
    d_dx_dh = -c * x / (rho + x)
    # ∂(dx/dt)/∂x
    d_dx_dx = mu * (1 - 2 * x / k) - c * h * rho / (rho + x) ** 2

    J = np.array([
        [d_dh, d_dx],
        [d_dx_dh, d_dx_dx]
    ])
    return J


if __name__ == "__main__":
    init_b = 1
    init_d = 1
    init_alpha = 0.001
    init_beta = 1
    init_mu = defaults.MU
    init_k = defaults.K
    init_c = defaults.C
    init_rho = defaults.RHO

    max_h = 10
    max_x = 8

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.38, right=0.95)

    # Initial vector field
    dh_dt = lambda h, x: dh_dt_with_params(h, x, init_b, init_rho, init_d, init_alpha, init_beta, init_k)
    dx_dt = lambda h, x: dx_dt_with_params(h, x, init_mu, init_k, init_c, init_rho)

    stream, sx_grid, sy_grid = vectors.plot_streamlines(dh_dt, dx_dt, ax, max_h, max_x)
    ax.set_xlabel('h')
    ax.set_ylabel('x')
    ax.set_title(f'Stream: mu={init_mu:.2f}, b={init_b:.2f}, d={init_d:.2f}, alpha={init_alpha:.2f}, beta={init_beta:.2f}, c={init_c:.2f}, k={init_k:.2f}, rho={init_rho:.2f}')
    
    # Find and plot fixed points
    params = (init_b, init_d, init_alpha, init_beta, init_k, init_mu, init_c, init_rho)
    # Generate a grid of guesses in the plotting region
    h_guesses = np.linspace(0.01, max_h, 5)
    x_guesses = np.linspace(0.01, max_x, 5)
    guesses = [(h, x) for h in h_guesses for x in x_guesses]
    fixed_points = find_fixed_points(params, guesses)
    plotted = set()
    for h_star, x_star in fixed_points:
        # Avoid plotting duplicates
        key = (round(h_star, 3), round(x_star, 3))
        if key in plotted:
            continue
        plotted.add(key)
        J = jacobian(h_star, x_star, params)
        eigs = np.linalg.eigvals(J)
        stable = np.all(np.real(eigs) < 0)
        color = 'green' if stable else 'red'
        ax.plot(h_star, x_star, 'o', color=color, markersize=10, markeredgecolor='black')

    # Slider axes
    ax_grid = gridspec.GridSpec(8, 1, left=0.07, bottom=0.05, right=0.3, top=0.95, wspace=0.4, hspace=2)
    ax_mu_slider = plt.subplot(ax_grid[0])
    ax_b_slider = plt.subplot(ax_grid[1])
    ax_d_slider = plt.subplot(ax_grid[2])
    ax_alpha_slider = plt.subplot(ax_grid[3])
    ax_beta_slider = plt.subplot(ax_grid[4])
    ax_c_slider = plt.subplot(ax_grid[5])
    ax_k_slider = plt.subplot(ax_grid[6])
    ax_rho_slider = plt.subplot(ax_grid[7])

    mu_slider = Slider(ax_mu_slider, 'mu (log)', -2, 2, valinit=log10(init_mu))
    b_slider = Slider(ax_b_slider, 'b (log)', -2, 2, valinit=log10(init_b))
    d_slider = Slider(ax_d_slider, 'd (log)', -2, 2, valinit=log10(init_b))
    alpha_slider = Slider(ax_alpha_slider, 'alpha (log)', -3, 0, valinit=log10(init_alpha))
    beta_slider = Slider(ax_beta_slider, 'beta (log)', 0, 3, valinit=log10(init_beta))
    c_slider = Slider(ax_c_slider, 'c (log)', -1, 1, valinit=log10(init_c))
    k_slider = Slider(ax_k_slider, 'k (log)', 0, 1, valinit=log10(init_k))
    rho_slider = Slider(ax_rho_slider, 'rho (log)', -1, 1, valinit=log10(init_rho))

    def update(val):
        mu = 10 ** mu_slider.val
        b = 10 ** b_slider.val
        d = 10 ** d_slider.val
        alpha = 10 ** alpha_slider.val
        beta = 10 ** beta_slider.val
        c = 10 ** c_slider.val
        k = 10 ** k_slider.val
        rho = 10 ** rho_slider.val

        ax.cla()
        new_dh_s = dh_dt_with_params(sx_grid, sy_grid, b, rho, d, alpha, beta, k)
        new_dx_s = dx_dt_with_params(sx_grid, sy_grid, mu, k, c, rho)
        ax.streamplot(sx_grid, sy_grid, new_dh_s, new_dx_s)
        ax.set_xlim(-0.1, max_h)
        ax.set_ylim(-0.1, max_x)
        ax.set_xlabel('h')
        ax.set_ylabel('x')
        ax.set_title(f'Stream: {mu=:.2f}, {b=:.2f}, {d=:.2f}, {alpha=:.2f}, {beta=:.2f}, {c=:.2f}, {k=:.2f}, {rho=:.2f}')

        # Find and plot fixed points for updated params
        params = (b, d, alpha, beta, k, mu, c, rho)
        h_guesses = np.linspace(0.01, max_h, 5)
        x_guesses = np.linspace(0.01, max_x, 5)
        guesses = [(h, x) for h in h_guesses for x in x_guesses]
        fixed_points = find_fixed_points(params, guesses)
        plotted = set()
        for h_star, x_star in fixed_points:
            key = (round(h_star, 3), round(x_star, 3))
            if key in plotted:
                continue
            plotted.add(key)
            J = jacobian(h_star, x_star, params)
            eigs = np.linalg.eigvals(J)
            stable = np.all(np.real(eigs) < 0)
            color = 'green' if stable else 'red'
            ax.plot(h_star, x_star, 'o', color=color, markersize=10, markeredgecolor='black')

        fig.canvas.draw_idle()

    mu_slider.on_changed(update)
    b_slider.on_changed(update)
    d_slider.on_changed(update)
    alpha_slider.on_changed(update)
    beta_slider.on_changed(update)
    c_slider.on_changed(update)
    k_slider.on_changed(update)
    rho_slider.on_changed(update)

    plt.show()
