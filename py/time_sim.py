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

#### FLOATING POINT SIM
flt_x = x
flt_cos = np.cos(x)

#### FIXED POINT SIM
fxp_x = quantize(x, 2)
fxp_cos = np.cos(fxp_x)

pce_rv_1 = np.random.uniform(low=-1, high=1, size=size)
pce_rv_2 = np.random.uniform(low=-1, high=1, size=size)
pce =    0.00260417 + 0.00520833*(1.5*np.square(pce_rv_2) - 0.5) + 0.125*(pce_rv_1 * pce_rv_2)
# pce = -0.00153255 + -0.0840605*(pce_rv_2) + -0.00306511*(1.5*np.square(pce_rv_2) - 0.5) + 0.00210546*(pce_rv_1) + -0.0735626*(pce_rv_1 + pce_rv_2)
noise_pwr = np.mean(np.square(flt_cos - fxp_cos))

print("flt_cos - fxp_cos: {0}".format(noise_pwr))
print("Predicted: {0}".format(np.mean(np.square(pce))))

plt.hist(flt_cos - fxp_cos, bins=1000, density=True)
plt.hist(pce, bins=1000, density=True)
plt.show()

exit()
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

size=100000

a = np.pi/2.0
x = a + np.random.uniform(low=-1.0, high=1.0, size=size)

y = np.sin(x)

pce_rv = np.random.uniform(low=-1, high=1, size=size)
# y_pce = 1.1 * pce_rv
y_pce = (5.0/6.0) - (1.0/3.0)*(1.5*np.square(pce_rv) - 0.5)


print("Actual var: {0}".format(np.mean(np.square(y))))
print("Predicted var: {0}".format(np.mean(np.square(y_pce))))


plt.hist(y, bins=1000, density=True)
plt.hist(y_pce, bins=1000, density=True)
plt.show()
exit()

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

iters = 15

x = np.sin(2.0 * np.pi * (1.0/15.0) * np.arange(0, iters))

state_y = np.zeros(iters)

a = np.zeros(iters)
y = np.zeros(iters)
y_d = np.zeros(iters)
y_d_c = np.zeros(iters)

fxp_x_c = np.zeros(iters)
fxp_y = np.zeros(iters)
fxp_y_d = np.zeros(iters)
fxp_y_d_c = np.zeros(iters)


decay = 0.9
c1 = decay
c2 = 1.0 - decay

for i in range(iters):
    a[i] = x[i] * c2
    
    if(i - 1 < 0):
        y_d[i] = 0
    else:
        y_d[i] = y[i-1]
    
    y_d_c[i] = y_d[i] * c1

    y[i] = a[i] + y_d_c[i]

plt.plot(x)
plt.plot(y)
plt.grid()
plt.tight_layout()
plt.show()