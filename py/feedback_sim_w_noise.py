import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

size = 100000
bins = 1000

iters = 50

u = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))
n1 = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))
n2 = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))
n3 = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))
n4 = np.random.uniform(low=-1.0, high=1.0, size=(iters, size))


state_u = np.random.uniform(low=-1.0, high=1.0, size=size)
state_y = np.zeros(size)

a = np.zeros((iters, size))
a_n = np.zeros((iters, size))
b = np.zeros((iters, size))
c = np.zeros((iters, size))
c_n = np.zeros((iters, size))
d = np.zeros((iters, size))
y = np.zeros((iters, size))
y_d = np.zeros((iters, size))
y_d_c = np.zeros((iters, size))
y_d_c_n = np.zeros((iters, size))
e = np.zeros((iters, size))
e_n = np.zeros((iters, size))



for i in range(iters):
    a[i] = u[i] * 0.05
    a_n[i] = a[i] + n1[i]

    if (i - 1 < 0):
        b[i] = state_u
        y_d[i] = state_y
    else:
        b[i] = u[i-1]
        y_d[i] = y[i-1]
   
    e[i] = a_n[i] * y_d[i]
    e_n[i] = e[i] + n4[i]

    c_n[i] = c[i] + n2[i]

    y_d_c_n[i] = y_d_c[i] + n3[i]

    d[i] = c_n[i] + e_n[i]

    y[i] = d[i] + y_d_c_n[i]


print(np.var(y[iters-1]))
plt.hist(y[iters-1], bins=bins, density=True, label="Monte Carlo Sim")
plt.show()


    