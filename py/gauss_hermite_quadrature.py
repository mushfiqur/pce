import numpy as np
from scipy.special import eval_hermitenorm
import matplotlib.pyplot as plt

k = 100
x, w = np.polynomial.hermite.hermgauss(k)

def f(x):
    # return np.cos(x)
    return x

def expect_std(n):
    if n%2!=0:
        return 0
    else:
        return np.power(-1.0, n/2.0)*np.exp(-0.5)/np.math.factorial(n)

sum = 0.0

idx = 1
scale =  np.math.factorial(idx) * np.sqrt(np.pi)

mu = 0.5
var = 4.0
sigma = np.sqrt(var)

for i in range(len(w)):
    point = np.sqrt(2)*sigma*x[i] + mu
    sum += w[i] * f(point) * eval_hermitenorm(idx, np.sqrt(2)*x[i])

sum /= scale

# print("expect c{0}: {1}".format(idx, expect_std(idx)))
print("       c{0}: {1}".format(idx, sum))

actual = np.random.normal(loc=mu, scale=sigma, size=100000)
predict = sum * np.random.normal(loc=0.0, scale=1.0, size=100000) + mu

plt.hist(actual, bins=1000, density=True)
plt.hist(predict, bins=1000, density=True)
plt.show()