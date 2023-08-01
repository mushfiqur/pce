import matplotlib.pyplot as plt
import numpy as np

def quantize(f_val, bits):
    # return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

size = 100000
bitwidth = 48

k = 1
N = 15
K_p = 0.2667
K_i = 0.0178
K_0 = 1

input_signal = np.zeros(size)
fxp_input_signal = np.zeros(size)

integrator_out = 0
phase_estimate = np.zeros(size)
e_D = np.zeros(size)
e_F = np.zeros(size)
sin_out = np.zeros(size)
cos_out = np.ones(size)

fxp_integrator_out = 0
fxp_phase_estimate = np.zeros(size)
fxp_e_D = np.zeros(size)
fxp_e_F = np.zeros(size)
fxp_sin_out = np.zeros(size)
fxp_cos_out = np.ones(size)

#### FLOATING POINT SIM
# for n in range(99):
    # input_signal[n] = np.cos(2*np.pi*(k/N)*n + np.pi)
input_signal = np.random.uniform(low=-1, high=1, size=100000)

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
        phase_estimate[n] = phase_estimate[n-1] + K_0 * e_F[n]
    else:
        phase_estimate[n] = K_0 * e_F[n]

    sin_out[n] = -np.sin(2*np.pi*(k/N)*(n) + phase_estimate[n])
    cos_out[n] = np.cos(2*np.pi*(k/N)*(n) + phase_estimate[n])

#### FIXED POINT SIM
for n in range(99):
    fxp_input_signal[n] = quantize(input_signal[n], bitwidth)

for n in range(99):
    if(n - 1 > 0):
        fxp_e_D[n] = quantize((fxp_input_signal[n] * fxp_sin_out[n-1]), bitwidth)
    else:
        fxp_e_D[n] = 0.0

    #loop filter
    fxp_integrator_out += K_i * fxp_e_D[n]
    fxp_e_F[n] = K_p * fxp_e_D[n] + fxp_integrator_out

    #NCO
    if(n - 1 > 0):
        fxp_phase_estimate[n] = fxp_phase_estimate[n-1] + K_0 * fxp_e_F[n]
    else:
        fxp_phase_estimate[n] = K_0 * fxp_e_F[n]

    fxp_sin_out[n] = -np.sin(2*np.pi*(k/N)*(n) + fxp_phase_estimate[n])
    fxp_cos_out[n] = np.cos(2*np.pi*(k/N)*(n) + fxp_phase_estimate[n])


#### RESULTS
print("{0}".format(np.mean(np.square(e_F - fxp_e_F))))

plt.plot(e_F - fxp_e_F)
plt.show()
