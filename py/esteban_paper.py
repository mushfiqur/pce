'''
Example taken from: 
L. Esteban, J. Antonio López Martín and A. Regadío, "Round-off noise estimation of
 fixed-point algorithms using Modified Affine Arithmetic and Legendre Polynomials,"
 2020 XXXV Conference on Design of Circuits and Integrated Systems (DCIS), 
 Segovia, Spain, 2020, pp. 1-6, doi: 10.1109/DCIS51330.2020.9268668.

https://ieeexplore.ieee.org/document/9268668
'''

import numpy as np
import matplotlib.pyplot as plt
import time
np.random.seed(int(time.time()))

def quantize(f_val, bits):
    return ((f_val * np.power(2.0, bits))).round() / np.power(2.0, bits)
    # return (f_val * np.power(2.0, bits)).astype(int) / np.power(2.0, bits)

bitwidth = 10

size = 1000000

eps = np.random.uniform(low=-1, high=1, size=size)

#### FLOATING POINT SIM
flt_a = 10.0 * eps
flt_b = 10.0 * eps
flt_c = 1.0 * eps
flt_d = 1.0 * eps

flt_e = flt_a * flt_b
flt_f = flt_c * flt_d

flt_g = flt_e - flt_f

#### FIXED POINT SIM
fxp_a = 10.0*quantize(eps, 8)
fxp_b = 10.0*quantize(eps, 8)
fxp_c =  1.0*quantize(eps, 8)
fxp_d =  1.0*quantize(eps, 8)

fxp_e = quantize(fxp_a * fxp_b, 1)
fxp_f = quantize(fxp_c * fxp_d, 1)

fxp_g = quantize(fxp_e - fxp_f, 11)

#### RESULTS

print("flt_g - fxp_g: {0}".format(np.mean(np.square(flt_g - fxp_g))))
