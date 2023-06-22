import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

size = 1000000
bins = 1000

x = np.random.uniform(low=-1.0, high=1.0, size=size)
# sim_mc = x * x
# sim_mc = x * x * x
sim_mc = x * x * x * x

x = np.random.uniform(low=-1.0, high=1.0, size=size)
# sim_pce = (1.0/3.0) + (2.0/3.0)*((3.0/2.0)*np.square(x) - (1.0/2.0))
# sim_pce = 0.6 * x
sim_pce = 0.2 + 0.4*((3.0/2.0)*np.square(x) - (1.0/2.0))


plt.hist(sim_mc, bins=bins, density=True, label="Monte Carlo Sim")
plt.hist(sim_pce, bins=bins, density=True, label="PCE Sim")

mean_mc = np.mean(sim_mc)
mean_pce = np.mean(sim_pce)
mean_err = round(np.abs((mean_mc - mean_pce) / mean_mc) * 100.0, 5)

var_mc = np.var(sim_mc)
var_pce = np.var(sim_pce)
var_err = round(np.abs((var_mc - var_pce) / var_mc) * 100.0, 5)

print("Monte Carlo Mean: {0}\nEstimated Mean: {1}\nError: {2}%\n".format(mean_mc, mean_pce, mean_err))
print("Monte Carlo Var: {0}\nEstimated Var: {1}\nError: {2}%\n".format(var_mc, var_pce, var_err))

plt.show()

