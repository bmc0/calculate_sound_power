#!/usr/bin/env python

import sys, getopt, math, csv, array

# Equations from https://www.princeton.edu/3D3A/Publications/Tylka_3D3A_DICalculation.pdf

def omega(theta):
	return 2.0 * math.pi * (1.0 - math.cos(theta))

# Note: delta must be constant
def weight(delta, n, l):
	if n == 0 or n == l - 1:
		return 1.0 / (4.0 * math.pi) * omega(delta / 2.0)  # Note: equation in paper expected 0° and 180° to be included in both orbits
	else:
		return 1.0 / (4.0 * math.pi) * ((omega(n * delta + delta / 2.0) - omega(n * delta - delta / 2.0)) / 4.0)

def get_weight_sum(weights):
	w_sum = 0.0
	w_sum += weights[0] + weights[len(weights) - 1]
	for i in range(1, len(weights) - 1):
		w_sum += weights[i] * 4.0
	return w_sum

def read_data(filename, col, correction):
	print("Reading " + filename + "...")
	values = array.array('d')
	with open(filename, newline="") as f:
		r = csv.reader(f, delimiter=' ')
		for row in r:
			if row[0] != "*":
				values.append(float(row[col]))
	if correction != None:
		for i in range(0, len(correction)):
			values[i] += correction[i]
	return values

def add_partial_orbit(power_spectrum, correction, weights, start, end, interval, file_prefix, is_negative):
	for i in range(start, end):
		v = read_data("{0:s}{1:g}.txt".format(file_prefix, (-interval if is_negative else interval) * i), 1, correction)
		for j in range(0, len(power_spectrum)):
			power_spectrum[j] += math.pow(10.0, v[j] / 10.0) * weights[i]

def write_power_spectrum(filename, freqs, power_spectrum):
	print("Writing " + filename + "...")
	with open(filename, 'w', newline="") as f:
		w = csv.writer(f, delimiter=' ')
		for i in range(0, len(power_spectrum)):
			w.writerow([freqs[i], 10.0 * math.log10(power_spectrum[i]), 0.0])

default_interval     = 10.
default_di_offset_db = 50.

def usage():
	print("Usage: {0:s} [OPTION...]".format(sys.argv[0]))
	print("")
	print("Options:")
	print("    --interval=n         measurement interval in degrees (default: {0:g})".format(default_interval))
	print("    --di-offset=n        DI curve offset in dB (default: {0:g})".format(default_di_offset_db))
	print("    --mirror-horizontal  use only positive angles for horizontal orbit")
	print("    --mirror-vertical    use only positive angles for vertical orbit")
	print("    -h, --help           show this help")

interval = default_interval
di_offset_db = default_di_offset_db
mirror_horiz = False
mirror_vert = False

try:
	optlist, args = getopt.gnu_getopt(sys.argv[1:], "h", ["interval=", "di-offset=", "mirror-horizontal", "mirror-vertical", "help"])
except getopt.GetoptError as e:
	print("error: " + str(e))
	usage()
	sys.exit(2)
for o, a in optlist:
	if o in ("-h", "--help"):
		usage()
		sys.exit()
	elif o == "--interval":
		interval = float(a)
	elif o == "--di-offset":
		di_offset_db = float(a)
	elif o == "--mirror-horizontal":
		mirror_horiz = True
	elif o == "--mirror-vertical":
		mirror_vert = True
	else:
		assert False, "unhandled option"

n_meas = int(round(360.0 / interval))  # per orbit
total_n_meas = n_meas * 2 - 2          # actual total
d = math.radians(360.0 / n_meas)       # delta
l = int(round(n_meas / 2.0)) + 1       # per half-orbit

print("Parameters:")
print("  Measurement interval:   {0:g}°".format(360.0 / n_meas))
print("  Number of measurements: {0:d} ({1:d} per orbit)".format(total_n_meas, n_meas))
print("  Mirror horizontal:      {0:s}".format(str(mirror_horiz)))
print("  Mirror vertical:        {0:s}".format(str(mirror_vert)))
print("")

## Build weight array
weights = array.array('d')
for i in range(0, l):
	weights.append(abs(weight(d, i, l)))
#print(weights)
#print("weight sum =", get_weight_sum(weights))

## Build frequency array
freqs = read_data("h0.txt", 0, None)
print("Files have {0:d} data points".format(len(freqs)))

## Initialize power spectrum array
power_spectrum = array.array('d')
for i in range(0, len(freqs)):
	power_spectrum.append(0.0)

## Load measurement correction
correction = None
try:
	correction = read_data("correction.txt", 1, None)
except OSError:
	print("correction.txt not present")
else:
	print("Using correction.txt")

## Read horizontal orbit data
add_partial_orbit(power_spectrum, correction, weights, 0, l, interval, "h", False)
add_partial_orbit(power_spectrum, correction, weights, 1, l - 1, interval, "h", not mirror_horiz)

## Read vertical orbit data
add_partial_orbit(power_spectrum, correction, weights, 1, l - 1, interval, "v", False)
add_partial_orbit(power_spectrum, correction, weights, 1, l - 1, interval, "v", not mirror_vert)

write_power_spectrum("sound_power.txt", freqs, power_spectrum)

## Read listening_axis.txt and calculate directivity index
try:
	listening_axis = read_data("listening_axis.txt", 1, correction)
except OSError:
	print("Failed to open listening_axis.txt; DI will not be calculated")
else:
	for i in range(0, len(power_spectrum)):
		listening_axis[i] = math.pow(10.0, listening_axis[i] / 10.0)
	write_power_spectrum("corrected_listening_axis.txt", freqs, listening_axis)
	di_offset = math.pow(10.0, float(di_offset_db) / 10.0)
	for i in range(0, len(power_spectrum)):
		listening_axis[i] = listening_axis[i] / power_spectrum[i] * di_offset
	write_power_spectrum("directivity_index.txt", freqs, listening_axis)
