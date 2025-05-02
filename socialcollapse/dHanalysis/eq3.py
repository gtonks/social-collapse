from math import log10
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from .dx_dt import dx_dt
from . import vectors


def dh_dt_with_params(h, x, mu_h, c_h):
    return mu_h * h - np.divide(c_h * h, x, out=np.full_like(h, np.inf), where=x!=0)


if __name__ == "__main__":
    init_mu_h = 1
    init_c_h = 1
    max_h = 10
    max_x = 6

    fig, axs = plt.subplots(1, 2)
    plt.subplots_adjust(bottom=0.25)

    # Initial vector field
    dh_dt = lambda h, x: dh_dt_with_params(h, x, init_mu_h, init_c_h)

    ax1 = axs[0]
    quiv, x_grid, y_grid = vectors.plot_field(dh_dt, dx_dt, ax1, max_h, max_x)
    ax1.set_xlabel('h')
    ax1.set_ylabel('x')
    ax1.set_title(f'Field: mu_h={init_mu_h:.2f}, c_h={init_c_h:.2f}')

    ax2 = axs[1]
    stream, sx_grid, sy_grid = vectors.plot_streamlines(dh_dt, dx_dt, ax2, max_h, max_x)
    ax2.set_xlabel('h')
    ax2.set_ylabel('x')
    ax2.set_title(f'Stream: mu_h={init_mu_h:.2f}, c_h={init_c_h:.2f}')
    
    # Slider axes
    ax_mu_slider = plt.axes([0.2, 0.15, 0.65, 0.03])
    ax_c_slider = plt.axes([0.2, 0.08, 0.65, 0.03])

    mu_slider = Slider(ax_mu_slider, 'mu_h (log)', -1, 1, valinit=log10(init_mu_h))  # log10(0.1) to log10(10)
    c_slider = Slider(ax_c_slider, 'c_h (log)', -2, 2, valinit=log10(init_c_h))

    def update(val):
        mu_h = 10 ** mu_slider.val
        c_h = c_slider.val

        # Update vector field
        new_dh = dh_dt_with_params(x_grid, y_grid, mu_h, c_h)
        new_dx = dx_dt(x_grid, y_grid)
        quiv.set_UVC(new_dh, new_dx)

        ax2.cla()
        new_dh_s = dh_dt_with_params(sx_grid, sy_grid, mu_h, c_h)
        new_dx_s = dx_dt(sx_grid, sy_grid)
        ax2.streamplot(sx_grid, sy_grid, new_dh_s, new_dx_s)
        ax2.set_xlabel('h')
        ax2.set_ylabel('x')

        # Update titles
        ax1.set_title(f'Field: mu_h={mu_h:.2f}, c_h={c_h:.2f}')
        ax2.set_title(f'Stream: mu_h={mu_h:.2f}, c_h={c_h:.2f}')

        fig.canvas.draw_idle()

    mu_slider.on_changed(update)
    c_slider.on_changed(update)

    plt.show()
