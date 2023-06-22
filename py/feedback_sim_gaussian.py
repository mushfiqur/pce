import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

size = 100000
bins = 1000

iters = 30

rv = np.random.normal(0.0, 1.4, (iters, size))

u = 0.01*(1) + 0.02*(rv) + 0.01*(np.square(rv) - 1.0) + 0.01*(np.power(rv, 3.0) - 3.0*rv) + 0.005*(np.power(rv, 4.0) - 6.0*np.power(rv, 2.0) + 3.0)

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


# pce = 0.0260458 + 0.630176*x + 0.0520915*(1.5*np.square(x) -0.5)
# pce = 0.0260458 + 0.630176*x + 0.0520915*(1.5*np.square(x) -0.5) + 0.00378871*(2.5*np.power(x, 3) + -1.5*x)
# var_pce = np.var(pce)

arr = []
for i in range(len(y)):
    arr.append(np.var(y[i]))
    print(np.var(y[i]))
print()
# exit()
var_pce = 0.00280855
# var_pce = 0.000865909
var_mc = np.var(y[iters-1])
var_err = round(np.abs((var_pce - var_mc) / var_mc) * 100.0, 5)

print("Monte Carlo Var: {0:.6f}\nEstimated Var: {1:.6f}\nError: {2}%\n".format(var_mc, var_pce, var_err))

plt.plot(arr)
plt.axhline(y=var_pce, color='red')
plt.show()

# plt.hist(y[iters-1], bins=bins, density=True, label="Monte Carlo Sim")
# plt.show()


    