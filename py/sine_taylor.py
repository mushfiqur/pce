import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    # return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

bitwidth = 10

size = 1000000

x = np.random.uniform(low=-1, high=1, size=size)
c1 = 1.0/6.0

#### FLOATING POINT SIM
flt_x = x
flt_x_2 = flt_x * flt_x
flt_x_3 = flt_x_2 * flt_x
flt_x_frac = flt_x_3 * c1
flt_sin_x = flt_x - flt_x_frac

#### FIXED POINT SIM
fxp_x = quantize(x, bitwidth)
fxp_x_2 = quantize(fxp_x * fxp_x, bitwidth)
fxp_x_3 = quantize(fxp_x_2 * fxp_x, bitwidth)
fxp_x_frac = quantize(fxp_x_3 * c1, bitwidth)
fxp_sin_x = fxp_x - fxp_x_frac

#### RESULTS
print("flt_sin_x - fxp_sin_x: {0}".format(np.mean(np.square(flt_sin_x - fxp_sin_x))))

