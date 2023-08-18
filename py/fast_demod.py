import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000

iters = 10

I = np.random.uniform(low=-1.0, high=1.0, size=size)
Q = np.random.uniform(low=-1.0, high=1.0, size=size)

c1 = 0.5

I_d1 = np.zeros(size)
I_d2 = np.zeros(size)

Q_d1 = np.zeros(size)
Q_d2 = np.zeros(size)

I_bar = np.zeros(size)
Q_bar = np.zeros(size)

a = np.zeros(size)
b = np.zeros(size)

pre_scale = np.zeros(size)
out = np.zeros(size)

for n in range(size):
    if(n - 1 >= 0):
        I_d1 = I[n-1]
        Q_d1 = Q[n-1]
    else:
        I_d1 = 0.0
        Q_d1 = 0.0

    if(n - 2 >= 0):
        I_d2 = I[n-2]
        Q_d2 = Q[n-2]
    else:
        I_d2 = 0.0
        Q_d2 = 0.0
    
    I_bar[n] = I[n] * I_d2
    Q_bar[n] = Q[n] * Q_d2
    
    a[n] = I_d1 * Q_bar[n]
    b[n] = Q_d1 * I_bar[n]

    pre_scale[n] = a[n] - b[n]
    out[n] = c1 * pre_scale[n]

##-----------------------------------------------------
x1 = np.random.uniform(low=-1.0, high=1.0, size=size)
x2 = np.random.uniform(low=-1.0, high=1.0, size=size)
x3 = np.random.uniform(low=-1.0, high=1.0, size=size)
x4 = np.random.uniform(low=-1.0, high=1.0, size=size)
x5 = np.random.uniform(low=-1.0, high=1.0, size=size)
x6 = np.random.uniform(low=-1.0, high=1.0, size=size)

# pce_expr = 0.5*x3*x5 - 0.5*x2*x6 + 0.5*x2*x4 - 0.5*x1*x5
# pce_pwr = np.mean(np.square(pce_expr))

pce_pwr = 1.0/9.0
out_pwr = np.mean(np.square(out))

print("MC - PCE: {0}".format(out_pwr - pce_pwr))