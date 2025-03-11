import random
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool


sys_size = 500
init_indv = sys_size
time_rate = 1
max_time = 2000
n_h_vals = 50
n_realizations = 1000

h_vals = np.linspace(0, 350, n_h_vals)

final_indv = np.zeros((n_h_vals, n_realizations))

def run_sim(h):
    cur_indv = init_indv
    t = random.expovariate(time_rate)
    while cur_indv != 0 and t < max_time:
        logistic_birth = 1
        logistic_death = 1 * cur_indv / sys_size
        functional_resp = h / (1 + cur_indv)
        sum_rxns = logistic_birth + logistic_death + functional_resp

        rand = random.random() * sum_rxns
        if rand < logistic_birth:
            cur_indv += 1
        else:
            cur_indv -= 1

        t += random.expovariate(time_rate)

    return cur_indv


with Pool() as pool:
    for h_idx, h in enumerate(h_vals):
        results_cur_h = pool.map(run_sim, [h] * n_realizations)
        final_indv[h_idx, :] = results_cur_h

for q in (0.025, 0.5, 0.975):
    final_indv_quants = np.quantile(final_indv, q, 1)
    plt.scatter(h_vals, final_indv_quants, label=f"Quantile {q}")

plt.legend()
plt.show()
