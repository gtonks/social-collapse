import random
import numpy as np
import matplotlib.pyplot as plt


sys_size = 500
init_indv = sys_size
time_rate = 1
max_time = 2000
n_h_vals = 50
n_realizations = 1000

h_vals = np.linspace(0, 350, n_h_vals)

# indv_timelines = []
# event_timeses = []
final_indv = np.zeros((n_h_vals, n_realizations))

for h_idx, h in enumerate(h_vals):
    for i in range(n_realizations):
        cur_indv = init_indv
        # indv_timeline = [ cur_indv ]

        # event_times = [ 0 ]
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
            # indv_timeline.append(cur_indv)

            # event_times.append(t)
            t += random.expovariate(time_rate)

        # indv_timeline.append(cur_indv)
        # event_times.append(max_time)

        # indv_timelines.append(indv_timeline)
        # event_timeses.append(event_times)
        final_indv[h_idx, i] = cur_indv

# print(cur_indv)

# def type_i(t, y): return y * (1 - y/sys_size) - y / (1 + y)
# sol = solve_ivp(type_i, [0, max_time], [ 10 ])

# for i in range(n_realizations):
#     plt.plot(event_timeses[i], indv_timelines[i])
# # plt.plot(sol.t, sol.y[0])
# plt.show()

for q in (0.025, 0.5, 0.975):
    final_indv_quants = np.quantile(final_indv, q, 1)
    plt.scatter(h_vals, final_indv_quants, label=f"Quantile {q}")

# final_indv_meds = np.median(final_indv, 1)
# plt.scatter(h_vals, final_indv_meds)

# confidence = 0.95
# q = (1 - confidence) / 2
# lower = np.quantile(final_indv, q, 1)
# plt.scatter(h_vals, lower)

# upper = np.quantile(final_indv, 1-q, 1)
# plt.scatter(h_vals, upper)

plt.legend()
plt.show()
# print(f"Mean individuals at t={max_time}: {final_indv.mean()}")
# print(f"Standard deviation: {final_indv.std()}")
