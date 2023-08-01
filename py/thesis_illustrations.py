import matplotlib.pyplot as plt
import numpy as np

size = 100000
bins = 1000

set_frac_bits = 4
mean = np.power(2.0, -1*set_frac_bits) / -2.0
var = np.power(2.0, -2*set_frac_bits) / 12.0
min = mean - np.sqrt(3.0*var)
max = mean + np.sqrt(3.0*var)
print(mean)
print(var)
exit()
# x = np.random.random(size=size)
x = np.random.normal(loc=0.0, scale=1.0, size=size)
x_q = (x * np.power(2, set_frac_bits)).astype(int)
# x_q = (x * np.power(2, set_frac_bits)).round()

noise = x - (x_q / np.power(2, set_frac_bits))

plt.hist(noise, bins=bins, density=True)
# plt.axvline(x=min, color='red')
# plt.axvline(x=max, color='red')
plt.axhline(y=16, color='red')
plt.tight_layout()
plt.show()