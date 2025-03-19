import time
import numpy as np
import matplotlib.pyplot as plt

from sardinyes_h import run_series, run_parallel


def compare_series_parallel():
    _, axs = plt.subplots(1, 2)

    n_samples = 10
    hs = np.linspace(0, 5, n_samples)
    
    start = time.time()
    equilibria = run_series(hs)
    end=time.time()
    print(f'Total time: {end-start}')
    
    axs[0].plot(hs, equilibria)
    axs[0].set_xlim(hs.min(), hs.max())
    axs[0].set_title('Series Result')
    axs[0].set_xlabel('H')
    axs[0].set_ylabel('Resource Density at Equilibrium')

    start = time.time()
    equilibria = run_parallel(hs)
    end=time.time()
    print(f'Total time: {end-start}')
    
    axs[1].plot(hs, equilibria)
    axs[1].set_xlim(hs.min(), hs.max())
    axs[1].set_title('Parallel Result')
    axs[1].set_xlabel('H')
    axs[1].set_ylabel('Resource Density at Equilibrium')

    plt.show()


if __name__ == '__main__':
    compare_series_parallel()