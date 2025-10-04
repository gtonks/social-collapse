import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    data = np.load("final_search_over_b.npz")
    bs = data['bs']
    sigmas = data['sigmas']
    x_steady = data['x_steady']
    H_steady = data['H_steady']
    x_final = data['x_final']
    H_final = data['H_final']
    n_trials = x_final.shape[2]

    # Compute probability of resource availability (H > 0) over trials
    x_available = np.count_nonzero(x_final > 0.01, axis=-1)
    probs = x_available / n_trials

    # Create meshgrid for plotting
    X, Y = np.meshgrid(sigmas, bs)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, probs, cmap='viridis')
    ax.set_xlabel('sigma')
    ax.set_ylabel('b')
    ax.set_zlabel('Probability of resource availability')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    plt.show()