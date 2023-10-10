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
# x = np.cos(2.0 * np.pi * (1.0/15.0) * np.arange(0, size) + np.random.uniform(low=0.0, high=2.0*np.pi))

state_y = np.zeros(size)

a = np.zeros(size)
y = np.zeros(size)
y_d = np.zeros(size)
y_d_c = np.zeros(size)

fxp_x_c = np.zeros(size)
fxp_y = np.zeros(size)
fxp_y_d = np.zeros(size)
fxp_y_d_c = np.zeros(size)

c1 = 0.5
c2 = 0.5

for i in range(size):
    a[i] = x[i] * c2
    
    if(i - 1 < 0):
        y_d[i] = 0
    else:
        y_d[i] = y[i-1]
    
    y_d_c[i] = y_d[i] * c1

    y[i] = a[i] + y_d_c[i]

#### FIXED POINT SIM
x_bitwidth    = 4
x_c_bitwidth  = 6
y_dc_bitwidth = 7

fxp_x = quantize(x, x_bitwidth)

for i in range(size):
    fxp_x_c[i] = quantize(fxp_x[i] * c2, x_c_bitwidth)
    
    if(i - 1 < 0):
        fxp_y_d[i] = 0
    else:
        fxp_y_d[i] = fxp_y[i-1]
    
    fxp_y_d_c[i] = quantize(fxp_y_d[i] * c1, y_dc_bitwidth)

    fxp_y[i] = fxp_x_c[i] + fxp_y_d_c[i]


#### RESULTS

sig_pwr = np.mean(np.square(y))
noise_pwr = np.mean(np.square(y - fxp_y))

print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))

'''
Matlab suggests 45 bits total
{
    x: 14
    x_c: 15
    y_d_c: 16
}
To achieve a 90.78 dB SNR

I suggest 46 bits total
{
    x: 15
    x_c: 16
    y_d_c: 16
}
To achieve a 94.87 dB SNR
'''