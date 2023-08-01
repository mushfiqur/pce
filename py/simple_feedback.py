import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000

iters = 10

x = np.random.uniform(low=-1.0, high=1.0, size=size)

state_y = np.zeros(size)

a = np.zeros((iters, size))
y = np.zeros((iters, size))
y_d = np.zeros((iters, size))
y_d_c = np.zeros((iters, size))

fxp_a = np.zeros((iters, size))
fxp_y = np.zeros((iters, size))
fxp_y_d = np.zeros((iters, size))
fxp_y_d_c = np.zeros((iters, size))

c1 = 0.5

for i in range(iters):
    a[i] = x
    
    if(i - 1 < 0):
        y_d[i] = state_y
    else:
        y_d[i] = y[i-1]
    
    y_d_c[i] = y_d[i] * c1

    y[i] = a[i] + y_d_c[i]

#### FIXED POINT SIM
fxp_x = quantize(x, 18)
for i in range(iters):
    fxp_a[i] = fxp_x
    
    if(i - 1 < 0):
        fxp_y_d[i] = state_y
    else:
        fxp_y_d[i] = fxp_y[i-1]
    
    fxp_y_d_c[i] = quantize(fxp_y_d[i] * c1, 17)

    fxp_y[i] = fxp_a[i] + fxp_y_d_c[i]


#### RESULTS
sig_pwr = 1.0/3.0
noise_pwr = np.mean(np.square(y[iters-1] - fxp_y[iters-1]))

print("flt_out - fxp_out: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))
