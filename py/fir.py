import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

iters = 1
size = 100000

# b0 = 0.1
# b1 = 0.2
# b2 = 0.3
# b3 = 0.4
# b4 = 0.5

b0 = 0.7640
b1 = 0.6218
b2 = 0.2269
b3 = 0.7705
b4 = 0.8326

fir_in = np.random.uniform(low=-1, high=1, size=(iters, size))
out = np.zeros((iters, size))
fxp_out = np.zeros((iters, size))

for i in range(iters):
    out[i] = b0 * fir_in[i]
    
    if(i - 1 >= 0):
        out[i] += b1 * fir_in[i-1]
    
    if(i - 2 >= 0):
        out[i] += b2 * fir_in[i-2]

    if(i - 3 >= 0):
        out[i] += b3 * fir_in[i-3]
    
    if(i - 4 >= 0):
        out[i] += b4 * fir_in[i-4]

#### FIXED POINT SIM

bitwidth = 1

fxp_fir_in = quantize(fir_in, 1)

for i in range(iters):
    fxp_out[i] = quantize(b0 * fxp_fir_in[i], bitwidth)
    
    if(i - 1 >= 0):
        out[i] += quantize(b1 * fxp_fir_in[i-1], bitwidth)
    
    if(i - 2 >= 0):
        out[i] += quantize(b2 * fxp_fir_in[i-2], bitwidth)

    if(i - 3 >= 0):
        out[i] += quantize(b3 * fxp_fir_in[i-3], bitwidth)
    
    if(i - 4 >= 0):
        out[i] += quantize(b4 * fxp_fir_in[i-4], bitwidth)

print("out - fxp_out: {0}".format(np.mean(np.square(out - fxp_out))))