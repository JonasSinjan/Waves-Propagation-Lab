# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 02:29:52 2018

@author: lsw
"""

import numpy as np
from scipy.optimize import leastsq
import pylab as plt
import csv
#np.set_printoptions(threshold=np.nan)

N = 4800 # number of data points
oldt = np.linspace(0, 8*np.pi, N)
t = oldt[2300:3600]
w = 2*np.pi/120
#data = 3.0*np.sin(t+0.001) + 0.5 + np.random.randn(N) # create artificial data with noise
data_list=[]
with open('120firstset.csv', newline='') as csvfile:
    data_row = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in data_row:
        print(row[1])
        data_list.append(float(row[1]))
data = np.asarray(data_list[2300:3600])
#print(data)
guess_mean = np.mean(data)
guess_std = 2*np.std(data)/(2**0.5)
guess_phase = 0.8*np.pi

# we'll use this to plot our first estimate. This might already be good enough for you
data_first_guess = guess_std*np.sin(t+guess_phase) + guess_mean

# Define the function to optimize, in this case, we want to minimize the difference
# between the actual data and our "guessed" parameters
optimize_func = lambda x: x[0]*np.sin(t+x[1]) + x[2] - data
est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

# recreate the fitted curve using the optimized parameters
data_fit = est_std*np.sin(t+est_phase) + est_mean

plt.plot(data, '.')
plt.plot(data_fit, label='after fitting')
#plt.plot(data_first_guess, label='first guess')
plt.legend()
plt.show()
print(est_phase)
print(est_std)
print(est_mean)