import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import stats

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151

in_fname = "/home/mushf/pce/samples/samples0.raw"
raw_data = np.fromfile(in_fname, dtype='uint8')
iq_data = (np.float32(raw_data) - 128.0)/128.0

rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

i_ds = i_filt[::10000]
# q_ds = q_filt[::rf_decim]

plt.hist(i_filt, bins=1000, density=True)
plt.show()
plt.tight_layout()
exit()


# out = stats.beta.fit(i_ds)

alpha = 1.1816456349440776
beta = 1.110501453170722
loc = -0.366895450310948
scale = 0.7060839558959622

x = np.linspace(stats.beta.ppf(0.01, alpha, beta), stats.beta.ppf(0.99, alpha, beta), 100)
y = stats.beta.pdf(x, a=alpha, b=beta, loc=loc, scale=scale)

plt.plot(x, y)
plt.show()

exit()
print(out)
exit()
plt.hist(i_filt, bins=1000, density=True)
plt.show()
plt.tight_layout()
exit()
