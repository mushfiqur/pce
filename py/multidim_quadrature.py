import numpy as np
from scipy.special import eval_legendre, eval_hermitenorm
import matplotlib.pyplot as plt

def f(theta_1, theta_2):
    # return np.cos(theta_1)*theta_2
    return theta_1*theta_2

k=20
x_leg, w_leg = np.polynomial.legendre.leggauss(k)
x_herm, w_herm = np.polynomial.hermite.hermgauss(k)

# scale = (2*idx_1 + 1) / (2.0 * np.math.factorial(idx_2) * np.sqrt(np.pi))
# scale = ((2*idx_1 + 1) / (2.0)) * ((2*idx_2 + 1) / (2.0))


a = -1
b = 1

mu = 0.0
var = 1.0
sigma = np.sqrt(var)

max_order = 2

mat = np.zeros((max_order+1, max_order+1))

for idx_1 in range(max_order+1):
    for idx_2 in range(max_order+1):
        # if(idx_1 + idx_2 <= max_order):
        if(True):
            coeff = 0.0
            scale = (2*idx_1 + 1) / (2.0 * np.math.factorial(idx_2) * np.sqrt(np.pi))
        
            for i in range(len(x_leg)):
                for j in range(len(x_herm)):
                    point_leg = ((b-a)/2.0)*x_leg[i] + (a+b)/2.0
                    point_herm = np.sqrt(2)*sigma*x_herm[j] + mu

                    coeff += (w_leg[i]*w_herm[j]) * \
                        f(point_leg, point_herm) * \
                        (eval_legendre(idx_1, x_leg[i]) * eval_hermitenorm(idx_2, np.sqrt(2)*x_herm[j]))
            
            coeff *= scale
            if(np.abs(coeff) > 1.0e-14):
                print("       c{0}_{1}: {2}".format(idx_1, idx_2, coeff))
                mat[idx_1][idx_2] = coeff
            else:
                print("       c{0}_{1}: {2}".format(idx_1, idx_2, 0.0))
                mat[idx_1][idx_2] = 0.0

        


size=100000
rv1 = np.random.uniform(low=a, high=b, size=size)
rv2 = np.random.normal(loc=mu, scale=sigma, size=size)

rv_a = np.random.uniform(low=-1, high=1, size=size)
rv_b = np.random.normal(loc=0.0, scale=1.0, size=size)

actual = f(rv1, rv2)

predict = np.zeros(size)

for i in range(max_order + 1):
    for j in range(max_order + 1):
        predict += mat[i][j]*eval_legendre(i, rv_a)*eval_hermitenorm(j, rv_b)

print("\nActual var: {0}".format(np.mean(np.square(actual))))
print("Predict var: {0}".format(np.mean(np.square(predict))))
error = ((np.mean(np.square(actual)) - np.mean(np.square(predict))) / np.mean(np.square(actual))) * 100
print("Error: {0:.2f}%".format(error))

# plt.hist(x=actual, bins=1000, density=True, alpha=1, label='actual')
# plt.hist(x=predict, bins=1000, density=True, alpha=0.8, label='pce')
# plt.legend()
# plt.tight_layout()
# plt.show()

exit()
