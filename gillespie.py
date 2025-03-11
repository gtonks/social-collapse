import random
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# TODO: Fix all parameters except H (ex. Figure 16.2), make fig like varying epsilon in Sardinyes


sys_size = 50
init_indv = 10
time_rate = 1
max_time = 100
n_realizations = 5

indv_timelines = []
event_timeses = []

for _ in range(n_realizations):
    cur_indv = init_indv
    indv_timeline = [ cur_indv ]

    event_times = [ 0 ]
    t = random.expovariate(time_rate)
    while cur_indv != 0 and t < max_time:
        logistic_birth = 1
        logistic_death = cur_indv / sys_size
        functional_resp = (1 + cur_indv)**-1
        sum_rxns = logistic_birth + logistic_death + functional_resp

        rand = random.random() * sum_rxns
        if rand < logistic_birth:
            cur_indv += 1
        else:
            cur_indv -= 1
        indv_timeline.append(cur_indv)

        event_times.append(t)
        t += random.expovariate(time_rate)

    indv_timeline.append(cur_indv)
    event_times.append(max_time)

    indv_timelines.append(indv_timeline)
    event_timeses.append(event_times)

# print(cur_indv)

D = 2
B = (sys_size + D - 1) / (sys_size * sys_size + D * sys_size)
def type_i(t, y): return y - B * y * y - y / (D + y)
sol = solve_ivp(type_i, [0, max_time], [ 10 ])

for i in range(n_realizations):
    plt.plot(event_timeses[i], indv_timelines[i])
plt.plot(sol.t, sol.y[0])
plt.show()
