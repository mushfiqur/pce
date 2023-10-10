import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy import stats
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000
test_coeffs = np.loadtxt("/home/mushf/pce/filters/rf_filt.txt")
# test_coeffs = signal.firwin(51, 100e3/(2.4e6/2), window=('hann'))
sum_coeffs_sqr = np.sum(np.square(test_coeffs))

rv1 = np.random.uniform(low=-1, high=1, size=size)
# input = 2.3 + 3.4*rv1 + 5.2*(1.5*rv1*rv1 - 0.5)
input = 1.0*rv1

flt_out = np.zeros(size)

# same as signal.lfilter(test_coeffs, 1.0, input )
for curr_timestamp in range(size):
    for k in range(len(test_coeffs)):
        if(curr_timestamp - k >= 0):
            flt_out[curr_timestamp] += test_coeffs[k] * input[curr_timestamp - k]


# mean = np.sum(test_coeffs) * np.mean(input)
# var = sum_coeffs_sqr * np.mean(np.square(input - np.mean(input)))
# points = np.linspace(stats.norm.ppf(0.001,loc=mean,scale=np.sqrt(var)), stats.norm.ppf(0.9999,loc=mean,scale=np.sqrt(var)),100)
# pdf = stats.norm.pdf(points,loc=mean,scale=np.sqrt(var))
# plt.hist(flt_out, bins=1000, density=True, alpha=1.0, label="actual")
# plt.plot(points, pdf, color='r', label="predicted")
# plt.tight_layout()
# plt.legend()
# plt.show()
# exit()


#### FIXED POINT SIM
fxp_input = quantize(input, 19)
fxp_out = np.zeros(size)

print("Doing fixed point sim...")
for curr_timestamp in range(size):
    for k in range(len(test_coeffs)):
        if(curr_timestamp - k >= 0):
            fxp_out[curr_timestamp] += quantize(test_coeffs[k] * fxp_input[curr_timestamp - k], 21)
            # fxp_out[curr_timestamp] += test_coeffs[k] * fxp_input[curr_timestamp - k]


#### RESULTS
sig_pwr = np.mean(np.square(flt_out))
noise_pwr = np.mean(np.square(flt_out - fxp_out))

print("Sig Pwr: {0}".format(sig_pwr))
# print("Noise pwr: {0}".format(noise_pwr))
print("flt_sin_x - fxp_sin_x: {0}".format(noise_pwr))
# print("                       {0}".format(len(test_coeffs) * (1.0/12.0)*(np.power(2.0, -2.0*bitwidth))))
print("SNR_Actual (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))
# print("SNR_Predicted (dB): {0}".format(10.0 * np.log10(sig_pwr / (len(test_coeffs) * (1.0/12.0)*(np.power(2.0, -2.0*bitwidth))))))


exit()

mean = np.sum(test_coeffs) * np.mean(input)
var = sum_coeffs_sqr * np.mean(np.square(input - np.mean(input)))

dist = np.random.normal(loc=mean, scale=np.sqrt(var), size=100000)


# PLOTTING
points = np.linspace(stats.norm.ppf(0.001,loc=mean,scale=np.sqrt(var)), stats.norm.ppf(0.9999,loc=mean,scale=np.sqrt(var)),100)
pdf = stats.norm.pdf(points,loc=mean,scale=np.sqrt(var))
plt.hist(flt_out, bins=1000, density=True, alpha=1.0, label="actual")
plt.plot(points, pdf, color='r', label="predicted")
plt.tight_layout()
plt.legend()
plt.show()
