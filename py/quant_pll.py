import matplotlib.pyplot as plt
import numpy as np
np.seterr(all='raise')

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100

k = 5
N = 15
K_p = 0.2667
K_i = 0.0178
K_0 = 1

input_signal = np.zeros(size)
fxp_input_signal = np.zeros(size)

integrator_out = 0.0
phase_estimate = np.zeros(size)
e_D = np.zeros(size)
e_F = np.zeros(size)
sin_out = np.zeros(size)
cos_out = np.ones(size)

fxp_integrator_out = 0.0
fxp_phase_estimate = np.zeros(size)
fxp_e_D = np.zeros(size)
fxp_e_F = np.zeros(size)
fxp_sin_out = np.zeros(size)
fxp_cos_out = np.ones(size)

#### FLOATING POINT SIM
# input_signal = np.sin(2*np.pi*(k/N)*np.arange(0, size)) + np.random.uniform(low=-1, high=1, size=size)
snr = 5.0
end = np.sqrt(3.0 * (0.5 / np.power(10.0, snr/10.0)))
# input_signal = np.sin(2*np.pi*(k/N)*np.arange(0, size)) + np.random.uniform(low=-end, high=end, size = size)
input_signal = np.sin(2*np.pi*(k/N)*np.arange(0, size))

for n in range(len(input_signal)):
    if(n - 1 > 0):
        e_D[n] = (input_signal[n] * sin_out[n-1])
    else:
        e_D[n] = 0.0

    #loop filter
    integrator_out += K_i * e_D[n]
    e_F[n] = K_p * e_D[n] + integrator_out

    #NCO
    if(n - 1 > 0):
        phase_estimate[n] = phase_estimate[n-1] + (K_0 * e_F[n]) + 2.0*np.pi*(k/N)
    else:
        phase_estimate[n] = K_0 * e_F[n] + 2.0*np.pi*(k/N)

    # sin_out[n] = -1.0*(x - (1.0/6.0)*np.power(x, 3))
    sin_out[n] = -1*np.sin(phase_estimate[n])
    cos_out[n] = np.cos(phase_estimate[n])

# plt.plot(input_signal)
# plt.plot(cos_out)
# plt.show()

# print(np.max(trig_arg))
# print("Actual pwr: {0}".format(np.mean(np.square(sin_out))))
# print("Predicted pwr: {0}".format(0.2))

#### FIXED POINT SIM
bitwidth = 1
for n in range(size):
    fxp_input_signal[n] = quantize(input_signal[n], 25)

for n in range(size):
    if(n - 1 > 0):
        fxp_e_D[n] = quantize((fxp_input_signal[n] * fxp_sin_out[n-1]), 33)
    else:
        fxp_e_D[n] = 0.0

    #loop filter
    fxp_integrator_out += quantize(K_i * fxp_e_D[n], 39)
    fxp_e_F[n] = quantize(K_p * fxp_e_D[n], 26) + fxp_integrator_out

    #NCO
    if(n - 1 > 0):
        fxp_phase_estimate[n] = fxp_phase_estimate[n-1] + K_0 * fxp_e_F[n] + 2.0*np.pi*(k/N)
    else:
        fxp_phase_estimate[n] = K_0 * fxp_e_F[n] + 2.0*np.pi*(k/N)

    fxp_sin_out[n] = -np.sin(fxp_phase_estimate[n])
    fxp_cos_out[n] =  np.cos(fxp_phase_estimate[n])


#### RESULTS
sig_pwr = np.mean(np.square(cos_out))
# sig_pwr = 0.5
noise_pwr = np.mean(np.square(cos_out - fxp_cos_out))

print("flt_out - fxp_out: {0}".format(noise_pwr))
print("SNR (dB): {0}".format(10.0 * np.log10(sig_pwr / noise_pwr)))

# fig, ax = plt.subplots(2)
# ax[0].plot(input_signal)
# ax[1].plot(cos_out)
# plt.tight_layout()
# plt.show()
