from collections.abc import Callable
import numpy as np
from matplotlib.axes import Axes


def plot_field(u: Callable[[np.ndarray, np.ndarray], np.ndarray], v: Callable[[np.ndarray, np.ndarray], np.ndarray], ax: Axes):
    step = 0.25
    max_x = 5
    max_y = 6
    x_range = np.arange(0, max_x+step, step)
    y_range = np.arange(0, max_y+step, step)
    x_grid, y_grid = np.meshgrid(x_range, y_range)
    us = u(x_grid, y_grid)
    vs = v(x_grid, y_grid)
    ax.quiver(x_grid, y_grid, us, vs)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)

def plot_direction(u: Callable[[np.ndarray, np.ndarray], np.ndarray], v: Callable[[np.ndarray, np.ndarray], np.ndarray], ax: Axes):
    step = 0.1
    max_x = 5
    max_y = 6
    x_range = np.arange(0, max_x+step, step)
    y_range = np.arange(0, max_y+step, step)
    x_grid, y_grid = np.meshgrid(x_range, y_range)
    us = u(x_grid, y_grid)
    vs = v(x_grid, y_grid)
    magnitudes = np.sqrt(us * us + vs * vs)
    u_dirs = np.divide(us, magnitudes, out=np.zeros_like(us), where=magnitudes!=0)
    i,j=np.where((0.09<x_grid) & (x_grid<0.11) & (2.89<y_grid) & (y_grid<2.91))
    print(u_dirs[i,j])
    v_dirs = np.divide(vs, magnitudes, out=np.zeros_like(vs), where=magnitudes!=0)
    ax.quiver(x_grid, y_grid, u_dirs, v_dirs)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)

def plot_streamlines(u: Callable[[np.ndarray, np.ndarray], np.ndarray], v: Callable[[np.ndarray, np.ndarray], np.ndarray], ax: Axes):
    max_x = 5
    max_y = 6
    x_range = np.linspace(0, max_x, 1000)
    y_range = np.linspace(0, max_y, 1000)
    x_grid, y_grid = np.meshgrid(x_range, y_range)
    us = u(x_grid, y_grid)
    vs = v(x_grid, y_grid)
    ax.streamplot(x_grid, y_grid, us, vs)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)