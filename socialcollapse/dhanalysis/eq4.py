"""
A sink at low r, a node at high r. Both in between. Beta has the opposite effect.
"""

from math import log10
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from socialcollapse import defaults
from .dx_dt import dx_dt
from . import vectors


def gamma_with_params(x, alpha, beta, k):
    return beta - (beta - alpha) * x / k


def dh_dt_with_params(h, x, r, alpha, beta, k):
    return h * (r - gamma_with_params(x, alpha, beta, k) * h)


if __name__ == "__main__":
    init_r = 1
    alpha = 0
    init_beta = 1
    k = defaults.K

    max_h = 10
    max_x = 6

    fig, axs = plt.subplots(1, 2)
    plt.subplots_adjust(bottom=0.25)

    # Initial vector field
    dh_dt = lambda h, x: dh_dt_with_params(h, x, init_r, alpha, init_beta, k)

    ax1 = axs[0]
    quiv, x_grid, y_grid = vectors.plot_field(dh_dt, dx_dt, ax1, max_h, max_x)
    ax1.set_xlabel('h')
    ax1.set_ylabel('x')
    ax1.set_title(f'Field: r={init_r:.2f}, beta={init_beta:.2f}')

    ax2 = axs[1]
    stream, sx_grid, sy_grid = vectors.plot_streamlines(dh_dt, dx_dt, ax2, max_h, max_x)
    ax2.set_xlabel('h')
    ax2.set_ylabel('x')
    ax2.set_title(f'Stream: r={init_r:.2f}, beta={init_beta:.2f}')
    
    # Slider axes
    ax_mu_slider = plt.axes([0.2, 0.15, 0.65, 0.03])
    ax_c_slider = plt.axes([0.2, 0.08, 0.65, 0.03])

    r_slider = Slider(ax_mu_slider, 'r (log)', -1, 1, valinit=log10(init_r))  # log10(0.1) to log10(10)
    beta_slider = Slider(ax_c_slider, 'beta (log)', -1, 1, valinit=log10(init_beta))

    def update(val):
        r = 10 ** r_slider.val
        beta = 10 ** beta_slider.val

        # Update vector field
        new_dh = dh_dt_with_params(x_grid, y_grid, r, alpha, beta, k)
        new_dx = dx_dt(x_grid, y_grid)
        quiv.set_UVC(new_dh, new_dx)

        ax2.cla()
        new_dh_s = dh_dt_with_params(sx_grid, sy_grid, r, alpha, beta, k)
        new_dx_s = dx_dt(sx_grid, sy_grid)
        ax2.streamplot(sx_grid, sy_grid, new_dh_s, new_dx_s)
        ax2.set_xlabel('h')
        ax2.set_ylabel('x')

        # Update titles
        ax1.set_title(f'Field: r={r:.2f}, beta={beta:.2f}')
        ax2.set_title(f'Stream: r={r:.2f}, beta={beta:.2f}')

        fig.canvas.draw_idle()

    r_slider.on_changed(update)
    beta_slider.on_changed(update)

    plt.show()
