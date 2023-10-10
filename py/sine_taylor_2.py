import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 1000000

x = np.random.uniform(low=0.0, high=2.0*np.pi, size=size)
x_0 = 2.0

c0 = x_0
c1 = -1.0/6.0*np.cos(x_0)
c2 = -1.0/2.0*np.sin(x_0)
c3 = np.cos(x_0)
c4 = np.sin(x_0)

#### FLOATING POINT SIM
flt_x = x
flt_d = flt_x - x_0

flt_x_2 = flt_d * flt_d
flt_x_3 = flt_x_2 * flt_d
flt_top = flt_x_3 * c1
flt_mid = flt_x_2 * c2
flt_bot = flt_d * c3
flt_sin_x = flt_top + flt_mid + flt_bot + c4


#### FIXED POINT SIM
fxp_x = quantize(x, 5)
fxp_d = fxp_x - x_0

fxp_x_2 = quantize(fxp_d * fxp_d, 1)
fxp_x_3 = quantize(fxp_x_2 * fxp_d, 1)

fxp_top = quantize(fxp_x_3 * c1, 4)
fxp_mid = quantize(fxp_x_2 * c2, 2)
fxp_bot = quantize(fxp_d * c3, 4)

fxp_sin_x = fxp_top + fxp_mid + fxp_bot + c4

#### RESULTS
sig_pwr = np.mean(np.square(flt_sin_x))
noise_pwr = np.mean(np.square(flt_sin_x - fxp_sin_x))

print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))
