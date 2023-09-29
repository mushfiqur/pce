import numpy as np
from scipy.special import eval_legendre
import matplotlib.pyplot as plt

k = 100
x, w = np.polynomial.legendre.leggauss(k)

def f(x):
    # return np.cos(x)
    return x

idx = 1
scale = (2.0*idx + 1.0)/2.0

a = 0.0
b = 2.0*np.pi

sum = 0.0
for i in range(len(w)):
    point = ((b-a)/2.0)*x[i] + (a+b)/2.0
    sum += w[i] * f(point) * eval_legendre(idx, x[i])

sum *= scale

# print("expect c{0}: {1}".format(idx, expect_std(idx)))
print("       c{0}: {1}".format(idx, sum))

actual = np.random.uniform(low=a, high=b, size=100000)
predict = sum * np.random.uniform(low=-1.0, high=1.0, size=100000) + (a+b)/2.0

plt.hist(actual, bins=1000, density=True)
plt.hist(predict, bins=1000, density=True)
plt.show()
