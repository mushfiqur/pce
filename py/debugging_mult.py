import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

'''
size = 1000000
x = np.random.uniform(low=-1, high=1, size=size)

c1 = 0.5

orig_var = np.var(x)
print("Original var: {0}".format(orig_var))
y = c1 * x
new_var = np.var(y)
print("New var: {0}".format(new_var))

print("Expected New var: {0}".format(c1 * orig_var))
print("PCE Expected Var: {0}".format(c1 * (1.0/3.0)))
exit()
'''
size = 100000
iters = 40

c1 = 0.5
out = np.zeros((iters, size))
out_d = np.zeros(size)

for i in range(iters):
    b = np.random.uniform(low=-1, high=1, size=size)
    if(i - 1 < 0):
        out_d = np.zeros(size)
    else:
        out_d = out[i-1]
    
    out_dc = c1 * out_d
    out[i] = b + out_dc

mc_var_arr = np.zeros(iters)
pce_var_arr = np.zeros(iters)

for i in range(len(out)):
    mc_var_arr[i] = np.var(out[i])
    pce_var_arr[i] = np.square((np.power(c1, i+1) - 1) / (c1 - 1)) * (1.0/3.0)

var_pce = np.square((np.power(c1, iters+1) - 1) / (c1 - 1)) * (1.0/3.0)
var_mc = np.var(out[iters-1])

print("Actual var: {0}".format(var_mc))
print("Estimated var: {0}".format(var_pce))

fig, ax = plt.subplots(ncols=2, nrows=1)
ax[0].plot(mc_var_arr)
ax[1].plot(pce_var_arr)

plt.plot(mc_var_arr, linestyle='dotted', color='red', label="Monte Carlo")
plt.plot(pce_var_arr, label="PCE")
plt.xlabel("Iteration")
plt.ylabel("Variance")
plt.tight_layout()
plt.legend()

plt.show()
exit()
'''
iters = 10000

b = 0.0
out = np.zeros(iters)
out_d = 0.0
out_dc = 0.0

c1 = 0.5

var_b = 0.0
mc_var_arr = np.zeros(iters)


for i in range(iters):
    b = np.random.uniform(low=-1, high=1)

    if(i - 1 < 0):
        out_d = 0.0
    else:
        out_d = out[i-1]
    
    out_dc = c1 * out_d
    out[i] = b + out_dc

## initially the PCE input is just x1
## Each additional iteration, the 1(x1) coeff is buffered and added to itself
## If the gain is 1.0, then each added iteration is is 1(x1) from the new input added to last state
## So after 5 iterations, it is 5(x1), for example
## Make all three models converge on this fact.
##      -> C++ Monte Carlo Sim
##      -> C++ PCE Sim
##      -> Python sims

pce_iters = 7
coeffs = np.zeros(pce_iters)

pwr = np.zeros(pce_iters)

for i in range(pce_iters):
    if(i - 1 < 0):
        coeffs[i] = 1.0
        pwr[i] = 1.0/3.0
    else:
        coeffs[i] = c1*coeffs[i-1]
        pwr[i] = pwr[i-1] + (np.square(coeffs[i]) * (1.0/3.0))
    

# plt.plot(pwr)

var_pce = np.sum(np.square(coeffs)) * (1.0/3.0)

print("MC: {0}".format(np.var(out)))
print("PCE: {0}".format(var_pce))
print("Conv: {0}".format(pwr[len(pwr) - 1]))

# plt.show()
exit()

'''

'''
P0: 1
P1: (x2)
P2: (x2)^(2) - 1
P3: (x1)
P4: (x1)(x2)
P5: (x1)^(2) - 1
'''

'''
iters = 2
size = 1000000

b = np.zeros(size)
out = np.zeros((iters,size))
out_d = np.zeros(size)
out_dc = np.zeros(size)

c1 = 1.5

for i in range(iters):
    b = np.random.uniform(low=-1, high=1, size=size)
    
    if(i - 1 < 0):
        out_d = np.zeros(size)
    else:
        out_d = out[i-1]

    out_dc = c1 * out_d
    out[i] = b + out_dc

## initially the PCE input is just x1
## Each additional iteration, the 1(x1) coeff is buffered and added to itself
## If the gain is 1.0, then each added iteration is is 1(x1) from the new input added to last state
## So after 5 iterations, it is 5(x1), for example
## Make all three models converge on this fact.
##      -> C++ Monte Carlo Sim
##      -> C++ PCE Sim
##      -> Python sims

coeffs = np.zeros(iters)

for i in range(iters):
    if(i - 1 < 0):
        coeffs[i] = 1.0
    else:
        coeffs[i] = c1*coeffs[i-1]

var_pce = np.sum(coeffs) * (1.0/3.0)

print("MC: {0}".format(np.var(out[iters-1])))
print("PCE: {0}".format(var_pce))

exit()
'''

## PCE gives an estimate of the mean and variances of a distribution, not the PDF itself
## Below example shows a LOT of correlation, yet it works. So IIR example should work

size = 100000
x = np.random.uniform(low=-1, high=1, size=size)

y_sim = np.random.uniform(low=-1, high=1, size=size) + \
    	np.random.uniform(low=-1, high=1, size=size) + \
        np.random.uniform(low=-1, high=1, size=size) + \
        np.random.uniform(low=-1, high=1, size=size) + \
        np.random.uniform(low=-1, high=1, size=size)

var_mc = np.var(y_sim)

var_pce = 5.0*(1.0/3.0)
print(var_mc)
print(var_pce)

# sqr = x * x						 ## All the nonlinearities of adding x over and over is represented by the x^2 term
# something = 25 * np.mean(sqr)


exit()

x1_0 = 0.87
x1_1 = 0.12
x1_2 = 0.06
x1_3 = 0.91
x1_4 = 0.48
x1_5 = 0.33

x2_0 = 0.53
x2_1 = 0.05
x2_2 = 0.80
x2_3 = 0.07
x2_4 = 0.60
x2_5 = 0.70

rv1 = np.random.normal(loc=0.0, scale=1.0, size=size)
rv2 = np.random.normal(loc=0.0, scale=1.0, size=size)

x1 = x1_0 + x1_1*(rv2) + x1_2*(np.power(rv2, 2.0) - 1) + x1_3*(rv1) + x1_4*(rv1*rv2) + x1_5*(np.power(rv1, 2.0) - 1)
x2 = x2_0 + x2_1*(rv2) + x2_2*(np.power(rv2, 2.0) - 1) + x2_3*(rv1) + x2_4*(rv1*rv2) + x2_5*(np.power(rv1, 2.0) - 1)

y_0 = 1.3768
y_1 = 0.8847
y_2 = 1.2138
y_3 = 1.9594
y_4 = 2.7383
y_5 = 2.0596

y_sim = x1 * x2

y_pce = y_0 + y_1*(rv2) + y_2*(np.power(rv2, 2.0) - 1) + y_3*(rv1) + y_4*(rv1*rv2) + y_5*(np.power(rv1, 2.0) - 1)

var_mc = np.var(y_sim)
var_pce = np.var(y_pce)
var_err = round(np.abs((var_mc - var_pce) / var_mc) * 100.0, 2)
print("Monte Carlo Var: {0}\nEstimated Var: {1}\nError: {2}%\n".format(var_mc, var_pce, var_err))
