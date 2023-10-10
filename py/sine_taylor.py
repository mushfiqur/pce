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
c1 = -1.0/6.0

#### FLOATING POINT SIM
flt_x = x
flt_x_2 = flt_x * flt_x
flt_x_3 = flt_x_2 * flt_x
flt_x_frac = flt_x_3 * c1
flt_sin_x = flt_x - flt_x_frac

# np.divide(flt_sin_x, flt_x, out=np.zeros_like(flt_sin_x), where=flt_x!=0)

#### FIXED POINT SIM
fxp_x = quantize(x, 7)                             ## matlab suggested: 15 (14) vs. my suggestion: 16 (15)
fxp_x_2 = quantize(fxp_x * fxp_x, 4)               ## matlab suggested: 16 (15) vs. my suggestion: 13 (12)
fxp_x_3 = quantize(fxp_x_2 * fxp_x, 3)             ## matlab suggested: 15 (14) vs. my suggestion: 15 (14)
fxp_x_frac = quantize(fxp_x_3 * c1, 5)             ## matlab suggested: 17 (17) vs. my suggestion: 17 (15)
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
sig_pwr = np.mean(np.square(flt_sin_x))
noise_pwr = np.mean(np.square(flt_sin_x - fxp_sin_x))

print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))

exit()
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

''''''''''''''' IMPORTANT '''''''''''''''
###### TAYLOR SERIES AROUND POINT a #####
size=100000

a = np.pi/2.0
# a = 0.0
x = a + np.random.uniform(low=-1.0, high=1.0, size=size)

d = x - a

top_a = d * d
top_b = top_a * d
top_c = top_b * (1.0/6.0)*np.cos(a)
mid_a = top_a * (-1.0/2.0)*np.sin(a)
bot_a = d * np.cos(a)

y = np.sin(a) + top_c + mid_a + bot_a
# y = np.sin(x)

pce_rv = np.random.uniform(low=-1, high=1, size=size)
# y_pce = 1.1 * pce_rv
y_pce = (5.0/6.0) - (1.0/3.0)*(1.5*np.square(pce_rv) - 0.5)


print("Actual var: {0}".format(np.mean(np.square(y))))
print("Predicted var: {0}".format(np.mean(np.square(y_pce))))


# plt.hist(y, bins=1000, density=True)
# plt.hist(y_pce, bins=1000, density=True)
# plt.show()
exit()