import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000

c1 = 1.0

#### FLOATING POINT SIM
flt_I   = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Id1 = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Id2 = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Q   = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Qd1 = np.random.uniform(low=-1.0, high=1.0, size=size)
flt_Qd2 = np.random.uniform(low=-1.0, high=1.0, size=size)

flt_Ibar = flt_I - flt_Id2
flt_Qbar = flt_Q - flt_Qd2

flt_a = flt_Qbar * flt_Id1
flt_b = flt_Ibar * flt_Qd1

flt_prescale = flt_a - flt_b

flt_out = flt_prescale

#### FIXED POINT SIM
fxp_I   = quantize(flt_I,   16)
fxp_Id1 = quantize(flt_Id1, 16)
fxp_Id2 = quantize(flt_Id2, 18)
fxp_Q   = quantize(flt_Q,   16)
fxp_Qd1 = quantize(flt_Qd1, 16)
fxp_Qd2 = quantize(flt_Qd2, 17)

fxp_Ibar = fxp_I - fxp_Id2
fxp_Qbar = fxp_Q - fxp_Qd2

fxp_a = quantize(fxp_Qbar * fxp_Id1, 18)
fxp_b = quantize(fxp_Ibar * fxp_Qd1, 18)

fxp_prescale = fxp_a - fxp_b

fxp_out = fxp_prescale

#### RESULTS
sig_pwr = np.mean(np.square(flt_out))
noise_pwr = np.mean(np.square(flt_out - fxp_out))

print("Sig pwr: {0}".format(sig_pwr))
print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))
