import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

iters = 10
size = 100000

input_I = np.random.uniform(low=-1, high=1, size=(iters, size))
input_Q = np.random.uniform(low=-1, high=1, size=(iters, size))

c1 = 0.5

I_d_1 = np.zeros((iters, size))
I_d_2 = np.zeros((iters, size))
Q_d_1 = np.zeros((iters, size))
Q_d_2 = np.zeros((iters, size))
I_bar = np.zeros((iters, size))
Q_bar = np.zeros((iters, size))
a = np.zeros((iters, size))
b = np.zeros((iters, size))
pre_scale = np.zeros((iters, size))
out = np.zeros((iters, size))

for i in range(iters):
    # input_I = input_I[i]
    # input_Q = input_Q[i]

    if(i - 1 >= 0):
        I_d_1[i] = input_I[i-1]
        Q_d_1[i] = input_Q[i-1]
    if(i - 2 >= 0):
        I_d_2[i] = input_I[i-2]
        Q_d_2[i] = input_Q[i-2]
    
    I_bar[i] = input_I[i] - I_d_2[i]
    Q_bar[i] = input_Q[i] - Q_d_2[i]
    
    a[i] = I_d_1[i] * Q_bar[i]
    b[i] = Q_d_1[i] * I_bar[i]
    
    pre_scale[i] = a[i] - b[i]
    out[i] = pre_scale[i] * c1


print("out pwr: {0}".format(np.mean(np.square(out[iters-1]))))


plt.hist(out[1], bins=1000, density=True)
plt.hist(out[6], bins=1000, density=True)
plt.hist(out[iters-1], bins=1000, density=True)
plt.tight_layout()
plt.show()

