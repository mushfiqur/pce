import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

size = 100000
bins = 1000

iters = 15

x = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))

state_x = np.zeros(size)
state_y = np.zeros(size)

a = np.zeros((iters, size))
y = np.zeros((iters, size))
y_d = np.zeros((iters, size))
y_d_c = np.zeros((iters, size))

c1 = 0.5
c2 = 1.0 - c1

for i in range(iters):
    a[i] = x[i] * c1
    
    if(i - 1 < 0):
        y_d[i] = state_y
    else:
        y_d[i] = y[i-1]
    
    y_d_c[i] = y_d[i] * c2

    y[i] = a[i] + y_d_c[i]
    


# pce = 0.0260458 + 0.630176*x + 0.0520915*(1.5*np.square(x) -0.5) + 0.00378871*(2.5*np.power(x, 3) + -1.5*x)
var_pce = 0.3333333

y_pce = 1 * x
# for i in range(len(d)):
    # print(np.var(d[i]))

# exit()

var_mc = np.var(a[iters-1])
# var_pce = 0.127308
var_err = round(np.abs((var_mc - var_pce) / var_mc) * 100.0, 5)

print("Monte Carlo Var: {0}\nEstimated Var: {1}\nError: {2}%\n".format(var_mc, var_pce, var_err))

plt.hist(y[iters-1], bins=bins, density=True, label="Monte Carlo Sim")
plt.hist(y_pce[iters-1], bins=bins, density=True, label="PCE")
plt.show()