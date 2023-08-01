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

fxp_x_c = np.zeros((iters, size))
fxp_y = np.zeros((iters, size))
fxp_y_d = np.zeros((iters, size))
fxp_y_d_c = np.zeros((iters, size))

c1 = 0.5
c2 = 0.5
x[0] = 1.0
for i in range(iters):
    a[i] = x[0] * c2
    
    if(i - 1 < 0):
        y_d[i] = 0
    else:
        y_d[i] = y[i-1]
    
    y_d_c[i] = y_d[i] * c1

    y[i] = a[i] + y_d_c[i]

#### FIXED POINT SIM
fxp_x = quantize(x, 3)

for i in range(iters):
    fxp_x_c[i] = quantize(fxp_x[0] * c2, 4)
    
    if(i - 1 < 0):
        fxp_y_d[i] = 0
    else:
        fxp_y_d[i] = fxp_y[i-1]
    
    fxp_y_d_c[i] = quantize(fxp_y_d[i] * c1, 5)

    fxp_y[i] = fxp_x_c[i] + fxp_y_d_c[i]


#### RESULTS
sig_pwr = 1.0/3.0
noise_pwr = np.mean(np.square(y[iters - 1000:] - fxp_y[iters - 1000:]))

print("flt_out - fxp_out: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))
