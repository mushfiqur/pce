import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

size = 100000
bins = 1000

iters = 15

u = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))

state_u = np.zeros(size)
state_y = np.zeros(size)

a = np.zeros((iters, size))
b = np.zeros((iters, size))
c = np.zeros((iters, size))
d = np.zeros((iters, size))
y = np.zeros((iters, size))
y_d = np.zeros((iters, size))
y_d_c = np.zeros((iters, size))
e = np.zeros((iters, size))

c1 = 0.05
c2 = 0.25
c3 = 0.6
c4 = 0.3

for i in range(iters):
    a[i] = u[i] * c1

    if (i - 1 < 0):
        b[i] = state_u
        y_d[i] = state_y
    else:
        b[i] = u[i-1]
        y_d[i] = y[i-1]

    e[i] = a[i] * y_d[i]

    c[i] = b[i] * c2

    y_d_c[i] = y_d[i] * c3

    d[i] = c[i] + e[i]

    y[i] = d[i] + y_d_c[i]


x = np.random.uniform(low=-1.0, high=1.0, size=size)
# pce = 0.0260458 + 0.630176*x + 0.0520915*(1.5*np.square(x) -0.5)
# pce = 0.0260458 + 0.630176*x + 0.0520915*(1.5*np.square(x) -0.5) + 0.00378871*(2.5*np.power(x, 3) + -1.5*x)
# var_pce = np.var(pce)


for i in range(len(d)):
    print(np.var(d[i]))

exit()

var_mc = np.var(y[iters-1])
var_pce = 0.146078
var_err = round(np.abs((var_mc - var_pce) / var_mc) * 100.0, 5)

print("Monte Carlo Var: {0}\nEstimated Var: {1}\nError: {2}%\n".format(var_mc, var_pce, var_err))

# plt.hist(y[iters-1], bins=bins, density=True, label="Monte Carlo Sim")
# plt.show()


    