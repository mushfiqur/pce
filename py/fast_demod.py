import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000

iters = 10

c1 = 0.5

#### FLOATING POINT SIM
flt_I = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Q = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_I_d1 = np.zeros(size)
flt_I_d2 = np.zeros(size)
flt_Q_d1 = np.zeros(size)
flt_Q_d2 = np.zeros(size)
flt_I_bar = np.zeros(size)
flt_Q_bar = np.zeros(size)
flt_a = np.zeros(size)
flt_b = np.zeros(size)
flt_pre_scale = np.zeros(size)
flt_out = np.zeros(size)

for n in range(size):
    if(n - 1 >= 0):
        flt_I_d1[n] = flt_I[n-1]
        flt_Q_d1[n] = flt_Q[n-1]
    else:
        flt_I_d1[n] = 0.0
        flt_Q_d1[n] = 0.0

    if(n - 2 >= 0):
        flt_I_d2[n] = flt_I[n-2]
        flt_Q_d2[n] = flt_Q[n-2]
    else:
        flt_I_d2[n] = 0.0
        flt_Q_d2[n] = 0.0
    
    flt_I_bar[n] = flt_I[n] - flt_I_d2[n]
    flt_Q_bar[n] = flt_Q[n] - flt_Q_d2[n]
    
    flt_a[n] = flt_I_d1[n] * flt_Q_bar[n]
    flt_b[n] = flt_Q_d1[n] * flt_I_bar[n]

    flt_pre_scale[n] = flt_a[n] - flt_b[n]
    flt_out[n] = flt_pre_scale[n]

#### FIXED POINT SIM
fxp_I = np.zeros(size)
fxp_Q = np.zeros(size)
fxp_I_d1 = np.zeros(size)
fxp_I_d2 = np.zeros(size)
fxp_Q_d1 = np.zeros(size)
fxp_Q_d2 = np.zeros(size)
fxp_I_bar = np.zeros(size)
fxp_Q_bar = np.zeros(size)
fxp_a = np.zeros(size)
fxp_b = np.zeros(size)
fxp_pre_scale = np.zeros(size)
fxp_out = np.zeros(size)

i_bitwidth   = 18
q_bitwidth   = 18
a_bitwidth   = 16
b_bitwidth   = 16

fxp_I = quantize(flt_I, i_bitwidth)
fxp_Q = quantize(flt_Q, q_bitwidth)

# fxp_I = flt_I
# fxp_Q = flt_Q

for n in range(size):
    if(n - 1 >= 0):
        fxp_I_d1[n] = fxp_I[n-1]
        fxp_Q_d1[n] = fxp_Q[n-1]
    else:
        fxp_I_d1[n] = 0.0
        fxp_Q_d1[n] = 0.0

    if(n - 2 >= 0):
        fxp_I_d2[n] = fxp_I[n-2]
        fxp_Q_d2[n] = fxp_Q[n-2]
    else:
        fxp_I_d2[n] = 0.0
        fxp_Q_d2[n] = 0.0
    
    fxp_I_bar[n] = fxp_I[n] - fxp_I_d2[n]
    fxp_Q_bar[n] = fxp_Q[n] - fxp_Q_d2[n]
    
    fxp_a[n] = quantize(fxp_I_d1[n] * fxp_Q_bar[n], a_bitwidth)
    fxp_b[n] = quantize(fxp_Q_d1[n] * fxp_I_bar[n], b_bitwidth)

    # fxp_a[n] = fxp_I_d1[n] * fxp_Q_bar[n]
    # fxp_b[n] = fxp_Q_d1[n] * fxp_I_bar[n]


    fxp_pre_scale[n] = fxp_a[n] - fxp_b[n]
    fxp_out[n] = fxp_pre_scale[n]
    # fxp_out[n] = c1 * fxp_pre_scale[n]

sig_pwr = np.mean(np.square(flt_out))
noise_pwr = np.mean(np.square(flt_out - fxp_out))

print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))


