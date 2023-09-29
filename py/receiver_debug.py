import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import signal
np.random.seed(int(time.time()))

size = 100000

rf_taps = 101
rf_Fs = 2.4e6
rf_Fc = 100e3
if_taps = 101
if_Fc = 16e3
if_Fs = 240e3

def filter(h, x):
    y = np.zeros(len(x))

    for n in range(len(x)):
        for k in range(len(h)):
            if(n-k >= 0):
                y[n] += h[k]*x[n-k]
    
    return y

def filt_sample_by_sample(h, x, curr_timestamp):
    y = 0.0
    
    for k in range(len(h)):
        if(curr_timestamp-k >= 0 and curr_timestamp < len(x)):
            y += h[k]*x[curr_timestamp-k]
    
    return y


#### FRONT END
# print("Generating Input")
dist_a = -1.0
dist_b =  1.0
input_I = np.random.uniform(low=dist_a, high=dist_b, size=size)
input_Q = np.random.uniform(low=dist_a, high=dist_b, size=size)

# print("Reading RF coeffs")
rf_coeff = np.loadtxt("/home/mushf/pce/filters/rf_filt.txt")

print("Filter I")
filtered_I = filter(h=rf_coeff, x=input_I)

print("Filter Q")
filtered_Q = filter(h=rf_coeff, x=input_Q)



pce = 0.333664*np.random.uniform(low=-1, high=1, size=size) + \
        0.333668*np.random.uniform(low=-1, high=1, size=size) + \
        0.333668*np.random.uniform(low=-1, high=1, size=size)

print("PCE VAR: {0}".format(np.mean(np.square(pce))))
print("ACT VAR: {0}".format(np.mean(np.square(filtered_I))))

plt.hist(filtered_I, bins=1000, density=True)
plt.hist(pce, bins=1000, density=True)
plt.show()

exit()
#### DEMODULATOR
c1 = 0.5
I_d_1 = np.zeros(size)
I_d_2 = np.zeros(size)
Q_d_1 = np.zeros(size)
Q_d_2 = np.zeros(size)
I_bar = np.zeros(size)
Q_bar = np.zeros(size)
a = np.zeros( size )
b = np.zeros( size )
pre_scale = np.zeros( size )
demod_out = np.zeros( size )


filtered_I = 0.274279*np.random.uniform(low=-1, high=1, size=100000)
filtered_Q = 0.274279*np.random.uniform(low=-1, high=1, size=100000)

print("Demodulating")
for i in range(size):
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

rv1 = np.random.uniform(low=-1, high=1, size=size)
rv2 = np.random.uniform(low=-1, high=1, size=size)
rv3 = np.random.uniform(low=-1, high=1, size=size)
rv4 = np.random.uniform(low=-1, high=1, size=size)
rv5 = np.random.uniform(low=-1, high=1, size=size)
rv6 = np.random.uniform(low=-1, high=1, size=size)

c = 0.0376144
pce = c*rv3*rv5 + c*rv3*rv4 + c*rv2*rv6 + c*rv1*rv6

print("PCE VAR: {0}".format(np.mean(np.square(pce))))
print("ACT VAR: {0}".format(np.mean(np.square(demod_out))))

# print("\nfiltered_I Pwr: {0}".format(np.mean(np.square(filtered_I))))
# print("demod_out Pwr: {0}".format(np.mean(np.square(demod_out))))
exit()

#### RDS CHANNEL FILTER
rds_channel_filter = np.loadtxt("/home/mushf/pce/filters/rds_channel_filt.txt")

print("Apply RDS Channel Filter")
rds_channel = filter(h=rds_channel_filter, x=demod_out)

#### CHANNEL SQUARED
rds_channel_squared = np.zeros(size)
print("Squaring")
rds_channel_squared = np.square(rds_channel)

#### RDS CARRIER FILTER
rds_carrier_filter = np.loadtxt("/home/mushf/pce/filters/rds_carrier_filt.txt")

print("Filter RDS Carrier")
rds_carrier =filter(h=rds_carrier_filter, x=rds_channel_squared)

#### PLL
print("Apply PLL")
for i in range(size):
    rds_carrier[i] += np.sin(2.0 * np.pi * (114000.0/240000.0) * i)

integrator_out = 0.0
phase_estimate = np.zeros(size)
e_D = np.zeros(size)
e_F = np.zeros(size)
sin_out = np.zeros(size)
cos_out = np.zeros(size)
trig_arg = np.zeros(size)

K_p = 0.002666
K_i = 3.555e-6
K_0 = 1

for n in range(size):
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

    trig_arg[n] = 2*np.pi*(114000.0/240000.0)*(n) + phase_estimate[n]

    sin_out[n] = -np.sin(trig_arg[n])
    cos_out[n] = np.cos(trig_arg[n])


#### MIXER
print("Mixing to Baseband")
rds_mixed = rds_channel * cos_out
mixer_gain = 2.0
rds_mixed = mixer_gain * rds_mixed

#### RDS BASEBAND FILTER
rds_baseband = np.zeros(size)
rds_baseband_filter = np.loadtxt("/home/mushf/pce/filters/rds_baseband_filt.txt")

print("Baseband filter")
# rds_baseband = signal.lfilter(rds_baseband_filter, 1.0, rds_mixed)
rds_baseband = filter(h=rds_baseband_filter, x=rds_mixed)

#### RDS RRC FILTER
rds_rrc_filter = np.loadtxt("/home/mushf/pce/filters/rds_rrc_filt.txt")

print("RRC Filter")
rds_rrc = filter(h=rds_rrc_filter, x=rds_baseband)

# print("demod_out Pwr: {0}".format(np.mean(np.square(demod_out))))


# print("cos_out pwr: {0}".format(np.mean(np.square(cos_out))))


# arr = np.zeros(iters)
# for i in range(iters):
#     arr[i] = np.mean(np.square(cos_out[i]))
# plt.plot(arr)
# plt.axhline(y=np.mean(np.square(cos_out)), color='red')
# plt.show()