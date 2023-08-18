import matplotlib.pyplot as plt
import numpy as np

k = 1
N = 15
K_p = 0.2667
K_i = 0.0178
K_0 = 1

input_signal = np.zeros(100)

integrator_out = 0
phase_estimate = np.zeros(100)
e_D = np.zeros(100)
e_F = np.zeros(100)
sin_out = np.zeros(100)
cos_out = np.ones(100)

trig_arg = np.zeros(100)

for n in range(99):
    input_signal[n] = np.cos(2*np.pi*(k/N)*n + np.pi)

for n in range(99):
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

    trig_arg[n] = 2*np.pi*(k/N)*(n) + phase_estimate[n]

    sin_out[n] = -np.sin(trig_arg[n])
    cos_out[n] = np.cos(trig_arg[n])



# Create a Figure
fig = plt.figure()

# Set up Axes
ax1 = fig.add_subplot(211)
ax1.plot(cos_out, label='PLL Output')
plt.grid()
ax1.plot(input_signal, label='Input Signal')
plt.legend()
ax1.set_title('Waveforms')

# Show the plot
#plt.show()

ax2 = fig.add_subplot(212)
ax2.plot(e_F)
plt.grid()
ax2.set_title('Filtered Error')
plt.tight_layout()
plt.show()
