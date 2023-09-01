import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal

def fmDemodArctan(I, Q, prev_phase = 0.0):
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

rf_taps = 101
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_decim = 10

if_taps = 101
if_Fc = 16e3
if_Fs = 240e3

audio_Fs = 48e3
audio_decim = 5

rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
audio_coeff = signal.firwin(if_taps,if_Fc/(if_Fs/2), window = ('hann'))

# with open('./filters/rf_coeffs', 'w') as file:
# 	for i in rf_coeff:
# 		file.write(str(i) + "\n")

# with open('./filters/audio_coeffs', 'w') as file:
# 	for i in audio_coeff:
# 		file.write(str(i) + "\n")

in_fname = "./data/samples0.raw"
raw_data = np.fromfile(in_fname, dtype='uint8', count=100000)
iq_data = (np.float32(raw_data) - 128.0)/128.0

## Input PCE coeffs is X = {0.0, 0.373, -0.26} or {0.0, 0.424, 0.0}

# rv = np.random.uniform(low=-1, high=1, size=100000)
# # dist = 0.373*rv + -0.26*(1.5*np.square(rv) - 0.5)
# dist = 0.424*rv
# print("  E[X]: {0} vs. {1}".format(np.mean(iq_data[::2]), np.mean(dist)))
# print("E[X^2]: {0} vs. {1}".format(np.mean(np.square(iq_data[::2])), np.mean(np.square(dist))))
# print("E[X^3]: {0} vs. {1}".format(np.mean(np.power(iq_data[::2], 3.0)), np.mean(np.power(dist, 3.0))))


i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[::2])
q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

i_filt = i_filt[::rf_decim]
q_filt = q_filt[::rf_decim]

fm_demod, state_phase = fmDemodArctan(i_filt, q_filt)

audio_data = signal.lfilter(audio_coeff, 1.0, fm_demod)

audio_data = audio_data[::audio_decim]

out_fname = "./data/fmMonoBlock.wav"
wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))