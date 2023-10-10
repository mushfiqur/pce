import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import signal
np.random.seed(int(time.time()))

size = 100000
test_coeffs = np.loadtxt("/home/mushf/pce/filters/test_filt.txt")
sum_coeffs_sqr = 0.0
for i in test_coeffs:
    sum_coeffs_sqr += i*i

rv1 = np.random.uniform(low=-1, high=1, size=size)

input = 2.3 + 1.2*rv1 + 0.3*(1.5*rv1*rv1 - 0.5)

out = signal.lfilter(test_coeffs, 1.0, rv1 )

plt.hist(out, bins=1000, density=True)
plt.show()

exit()

rv3 = np.random.uniform(low=-1, high=1, size=size)
rv4 = np.random.uniform(low=-1, high=1, size=size)
rv5 = np.random.uniform(low=-1, high=1, size=size)
rv6 = np.random.uniform(low=-1, high=1, size=size)
rv7 = np.random.uniform(low=-1, high=1, size=size)
rv8 = np.random.uniform(low=-1, high=1, size=size)
rv9 = np.random.uniform(low=-1, high=1, size=size)

demod_pce = 0.0374885*rv5*rv6 + (-0.0374885)*rv4*rv6 + (-0.0374885)*rv3*rv8 + (0.0374885)*rv3*rv7

x_i = 0.0374886
c = x_i * np.sqrt(sum_coeffs_sqr)

# print(c)
# print(0.0374855*np.sqrt(sum_coeffs_sqr))
# exit()

one_term_filtered_actual = signal.lfilter(test_coeffs, 1.0, 0.0374885*rv5*rv6)

alpha = (0.0374885/3.0)*np.sqrt(3.0*sum_coeffs_sqr)
# one_term_filtered_predicted = alpha*np.random.uniform(low=-1, high=1, size=size)
one_term_filtered_predicted = 0.0374855*np.sqrt(sum_coeffs_sqr)*rv5*rv6

total_filtered_actual = signal.lfilter(test_coeffs, 1.0, demod_pce)
# total_filtered_predicted = alpha*np.random.uniform(low=-1, high=1, size=size) - alpha*np.random.uniform(low=-1, high=1, size=size) - alpha*np.random.uniform(low=-1, high=1, size=size) + alpha*np.random.uniform(low=-1, high=1, size=size)
total_filtered_predicted = 0.0374855*np.sqrt(sum_coeffs_sqr)*rv5*rv6 - 0.0374855*np.sqrt(sum_coeffs_sqr)*rv4*rv6 - 0.0374855*np.sqrt(sum_coeffs_sqr)*rv3*rv8 + 0.0374855*np.sqrt(sum_coeffs_sqr)*rv3*rv7

one_term_actual_var = np.mean(np.square(one_term_filtered_actual))
one_term_predicted_var = np.mean(np.square(one_term_filtered_predicted))

# print("Actual: {0}".format(one_term_actual_var))
# print("Predicted: {0}".format(one_term_predicted_var))

print("Actual: {0}".format(np.mean(np.square(total_filtered_actual))))
print("Predicted: {0}".format(np.mean(np.square(total_filtered_predicted))))

# exit()
plt.hist(one_term_filtered_actual, bins=1000, density=True, alpha=1.0, label='actual dist')
plt.hist(one_term_filtered_predicted, bins=1000, density=True, alpha=0.8, label="predicted dist")
plt.show()