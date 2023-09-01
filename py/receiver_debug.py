import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import signal
np.random.seed(int(time.time()))

iters = 150
size = 10000

rf_taps = 101
rf_Fs = 2.4e6
rf_Fc = 100e3
if_taps = 101
if_Fc = 16e3
if_Fs = 240e3

#### FRONT END
input_I = np.random.uniform(low=-1, high=1, size=(iters, size))
input_Q = np.random.uniform(low=-1, high=1, size=(iters, size))

rf_coeff = np.loadtxt("/home/mushf/pce/filters/rf_filt.txt")


filtered_I = np.zeros((iters,size))
filtered_Q = np.zeros((iters,size))

for i in range(iters):
    filtered_I[i] = signal.lfilter(rf_coeff, 1.0, input_I[i])
    filtered_Q[i] = signal.lfilter(rf_coeff, 1.0, input_Q[i])

#### DEMODULATOR
c1 = 0.5
I_d_1 = np.zeros((iters, size))
I_d_2 = np.zeros((iters, size))
Q_d_1 = np.zeros((iters, size))
Q_d_2 = np.zeros((iters, size))
I_bar = np.zeros((iters, size))
Q_bar = np.zeros((iters, size))
a = np.zeros((iters, size))
b = np.zeros((iters, size))
pre_scale = np.zeros((iters, size))
demod_out = np.zeros((iters, size))

for i in range(iters):
    if(i - 1 >= 0):
        I_d_1[i] = filtered_I[i-1]
        Q_d_1[i] = filtered_Q[i-1]
    if(i - 2 >= 0):
        I_d_2[i] = filtered_I[i-2]
        Q_d_2[i] = filtered_Q[i-2]
    
    I_bar[i] = filtered_I[i] - I_d_2[i]
    Q_bar[i] = filtered_Q[i] - Q_d_2[i]
    
    a[i] = I_d_1[i] * Q_bar[i]
    b[i] = Q_d_1[i] * I_bar[i]
    
    pre_scale[i] = a[i] - b[i]
    demod_out[i] = pre_scale[i] * c1

#### RDS CHANNEL FILTER
rds_channel = np.zeros((iters, size))
rds_channel_filter = np.loadtxt("/home/mushf/pce/filters/rds_channel_filt.txt")

for i in range(iters):
    rds_channel[i] = signal.lfilter(rds_channel_filter, 1.0, demod_out[i])

#### CHANNEL SQUARED
rds_channel_squared = np.zeros((iters, size))
for i in range(iters):
    rds_channel_squared[i] = rds_channel[i] * rds_channel[i]

#### RDS CARRIER FILTER
rds_carrier = np.zeros((iters, size))
rds_carrier_filter = np.loadtxt("/home/mushf/pce/filters/rds_carrier_filt.txt")

for i in range(iters):
    rds_carrier[i] = signal.lfilter(rds_carrier_filter, 1.0, rds_channel_squared[i])

#### PLL

for i in range(iters):
    rds_carrier[i] += np.sin(2.0 * np.pi * (1.0/15.0) * i)

integrator_out = np.zeros(size)
phase_estimate = np.zeros((iters, size))
e_D = np.zeros((iters, size))
e_F = np.zeros((iters, size))
sin_out = np.zeros((iters, size))
cos_out = np.zeros((iters, size))
trig_arg = np.zeros((iters, size))

K_p = 0.2667
K_i = 0.0178
K_0 = 1

for n in range(iters):
    if(n - 1 > 0):
        e_D[n] = (rds_carrier[n] * sin_out[n-1])
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

    trig_arg[n] = 2*np.pi*(1.0/15.0)*(n) + phase_estimate[n]

    sin_out[n] = -np.sin(trig_arg[n])
    cos_out[n] = np.cos(trig_arg[n])

print("cos_out pwr: {0}".format(np.mean(np.square(cos_out))))


# arr = np.zeros(iters)
# for i in range(iters):
#     arr[i] = np.mean(np.square(cos_out[i]))
# plt.plot(arr)
# plt.axhline(y=np.mean(np.square(cos_out)), color='red')
# plt.show()

#### MIXER
rds_mixed = rds_channel * cos_out

#### RDS BASEBAND FILTER
rds_baseband = np.zeros((iters, size))
rds_baseband_filter = np.loadtxt("/home/mushf/pce/filters/rds_baseband_filt.txt")

for i in range(iters):
    rds_baseband[i] = signal.lfilter(rds_baseband_filter, 1.0, rds_mixed[i])

#### RDS RRC FILTER
rds_rrc = np.zeros((iters, size))
rds_rrc_filter = np.loadtxt("/home/mushf/pce/filters/rds_rrc_filt.txt")

for i in range(iters):
    rds_rrc[i] = signal.lfilter(rds_rrc_filter, 1.0, rds_baseband[i])

