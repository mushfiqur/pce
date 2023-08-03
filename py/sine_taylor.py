import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

bitwidth = 1

size = 1000000

x = np.random.uniform(low=-1, high=1, size=size)
c1 = 1.0/6.0

#### FLOATING POINT SIM
flt_x = x
flt_x_2 = flt_x * flt_x
flt_x_3 = flt_x_2 * flt_x
flt_x_frac = flt_x_3 * c1
flt_sin_x = flt_x - flt_x_frac

# np.divide(flt_sin_x, flt_x, out=np.zeros_like(flt_sin_x), where=flt_x!=0)

#### FIXED POINT SIM
fxp_x = quantize(x, 15)                             ## matlab suggested: 15 (14) vs. my suggestion: 16 (15)
fxp_x_2 = quantize(fxp_x * fxp_x, 12)               ## matlab suggested: 16 (15) vs. my suggestion: 13 (12)
fxp_x_3 = quantize(fxp_x_2 * fxp_x, 14)             ## matlab suggested: 15 (14) vs. my suggestion: 15 (14)
fxp_x_frac = quantize(fxp_x_3 * c1, 15)             ## matlab suggested: 17 (17) vs. my suggestion: 17 (15)
# fxp_sin_x = quantize( fxp_x - fxp_x_frac, 9)      ## tot:              63 (60)                    61 (56)
fxp_sin_x = fxp_x - fxp_x_frac

'''
(0.000985925) { 
        eps: 3, 
        x^2: 3, 
        x^3: 3, 
        x_{frac}: 7, 
        sin_x: 6
}
'''
#### RESULTS
sig_pwr = 1.0/3.0
noise_pwr = np.mean(np.square(flt_sin_x - fxp_sin_x))

print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))

'''
With 1000000 samples, Matlab suggested this config:
{
    x: 14
    x_2: 15
    x_3: 14
    x_frac: 17 
} (tot: 60 bits)
For a SNR of 91.5 dB

With a target of 92 dB, I suggest:
{ 
    x: 15, 
    x_2: 12, 
    x_3: 14, 
    x_frac: 15
} (tot: 56 bits)
which hits 92.7 dB.
'''