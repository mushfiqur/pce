import numpy as np
import matplotlib.pyplot as plt

n = 100
size=100000
dist = np.zeros(size)

a = -2
b = 10.0
for i in range(n):
	dist += np.random.uniform(low=a, high=b, size=100000)

mean = n*(a + b)/2.0
var = n*np.square(b-a)/12.0

print("Expected Mean: {0}".format(mean))
print("Expected Var: {0}".format(var))
print("--------")

print("Mean: {0}".format(np.mean(dist)))
print("Var: {0}".format(np.mean(np.square(dist - np.mean(dist)))))
