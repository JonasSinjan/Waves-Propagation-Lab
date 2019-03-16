# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 10:48:05 2018

@author: Jonas
"""

import numpy as np
import csv
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as spo
import scipy.integrate as spi
import scipy.misc as spm

ylist = []
with open('480.csv', newline='') as csvfile:
    data_row = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in data_row:
        if float(row[1]) < 10:
            pass
        else:
            ylist.append(float(row[1]))
T = 480
w = (2 * np.pi) / T
w_error = ((2 * np.pi) / (T ** 2)) * 4
N = len(ylist)  # number of data points
# 4800
# ]480
oldt = np.linspace(0, N / 10, N)
t = oldt[:4801]


def sinfit(x, A, B, C):
    return A + B * np.sin((w * x) + C)


def sinfit_plus_linear(x, A, B, C, D):
    return A + B * np.sin(w * x + C) + D * x


xc = t
yv = ylist[:4801]

plt.figure(0)
# plt.plot(xc,yv,'.')
plt.ylabel('Temperature (\u2103)')
plt.xlabel('Time (seconds)')

popt, pcov = spo.curve_fit(sinfit_plus_linear, xc, yv, p0=[50, 0.2, 0.3, 1.0])
# print (popt)
perr = np.sqrt(np.diag(pcov))
# print(perr)

newyv = [x - popt[0] for x in yv]
# dataplot, = plt.plot(xc,newyv,'.',label="Recorded Temperature")
dataplot, = plt.plot(xc, yv, '.', label="Recorded Temperature")


# while popt[2]<0:
#     popt[2]+=2*np.pi
# if w == (2*np.pi)/360:
#     popt[2]+=np.pi #360 was done out of phase

def sinfit_plus_linear_derivative(x):
    return popt[0] + popt[1] * np.sin(w * x + popt[2]) + popt[3] * x


yfit = [sinfit_plus_linear(x, popt[0], popt[1], popt[2], popt[3]) for x in xc]
newyfit = [x - popt[0] for x in yfit]
fitplot, = plt.plot(xc, yfit, label="Sinusoidal Fit")
# fitplot, = plt.plot(xc,newyfit,label="Sinusoidal Fit")


transfactor = popt[1] / 63.7


def damp(trans):
    return (w * 8.22 ** 2) / (2 * (np.log(np.abs(trans))) ** 2)


D = damp(transfactor)


# print (D)


def dp(p):
    if p < 0: p += 2 * np.pi
    return (w * 8.22 ** 2) / (2 * (p) ** 2)


# while popt[2]<0:
#     popt[2]+=2*np.pi
#
# lmao = (np.pi-popt[2]-2*np.pi)/((2*np.pi))
# print(lmao)
#
# print(dp(lmao+2*np.pi))
# print(popt[2])
# print(dp(popt[2]))
# yerr=np.sqrt(pcov[0]**2+(np.sin(w*x+popt[2])**2)*pcov[1]**2+(((popt[1]*np.cos(w*x+popt[2]))**2)*pcov[2]**2)+(x**2)*pcov[3]**2)


xval = []
for x1 in range(len(yfit)):
    yval = yfit[x1]
    avg = popt[0]
    yerr = np.sqrt(perr[0] ** 2 + (np.sin(w * x1 + popt[2]) ** 2) * perr[1] ** 2 + (
                ((popt[1] * np.cos(w * x1 + popt[2])) ** 2) * perr[2] ** 2) + ((x1 ** 2) * perr[3] ** 2))
    # yerrmax = max(yerr)
    if -yerr < yval - avg < yerr:  # 0.05 for 240, 0.01 for 60 and 120
        dt = spm.derivative(sinfit_plus_linear_derivative, x1 / 10, dx=1e-6)
        if T == 360:
            if dt < 0:  # < for 360 otherwise >
                xval.append(x1)
            else:
                pass
        else:
            if dt > 0:  # < for 360 otherwise >
                xval.append(x1)
            else:
                pass

if T == 60:
    xphasepoints = [x / 10 for x in xval[5:26]]
    xphasepoints2 = [(x / 10) ** 2 for x in xval[5:26]]
    meanx = np.average(xphasepoints)
    varx = np.average(xphasepoints2) - meanx ** 2
    meanerror = np.sqrt(varx / len(xphasepoints))
    print(meanx, meanerror, varx)
    print(max(xphasepoints), min(xphasepoints))
    phaseerr = np.sqrt(meanx ** 2 * w_error ** 2 + w ** 2 * meanerror ** 2)
    print(meanx * w, phaseerr, meanerror * w)
    Dphase = (w * 8.22 ** 2) / (2 * (w * meanx) ** 2)
    # Derr = np.sqrt((((8.22**2)/(w*meanx**2))**2)*0.05**2+((8.22**2)/(2*(meanx**2)*w**2))*((7E-3)**2)+(8.22**2/(w*meanx**3))*meanerror**2)
    Derr = np.sqrt(
        (((8.22) / (w * (meanx ** 2))) ** 2) * (0.05) ** 2 + (((8.22) ** 2) / (2 * (meanx ** 2) * w ** 2)) ** 2 * (
                    (w_error) ** 2) + ((8.22) ** 2 / (w * meanx ** 3)) ** 2 * meanerror ** 2)
    errpercent = (Derr / Dphase) * 100
    print(Dphase, Derr, errpercent)

if T == 120:
    xphasepoints = [x / 10 for x in xval[0:6]]
    xphasepoints2 = [(x / 10) ** 2 for x in xval[0:6]]
    meanx = np.average(xphasepoints)
    varx = np.average(xphasepoints2) - meanx ** 2
    meanerror = np.sqrt(varx / len(xphasepoints))
    print(meanx, meanerror, varx)
    print(max(xphasepoints), min(xphasepoints))
    phaseerr = np.sqrt(meanx ** 2 * w_error ** 2 + w ** 2 * meanerror ** 2)
    print(meanx * w, phaseerr, meanerror * w)
    Dphase = (w * 8.22 ** 2) / (2 * (w * meanx) ** 2)
    Derr = np.sqrt(
        (((8.22) / (w * (meanx ** 2))) ** 2) * (0.05) ** 2 + (((8.22) ** 2) / (2 * (meanx ** 2) * w ** 2)) ** 2 * (
                    (w_error) ** 2) + ((8.22) ** 2 / (w * meanx ** 3)) ** 2 * meanerror ** 2)
    errpercent = (Derr / Dphase) * 100
    print(Dphase, Derr, errpercent)

if T == 240:
    xphasepoints = [x / 10 for x in xval[0:2]]
    print(xphasepoints)
    xphasepoints2 = [(x / 10) ** 2 for x in xval[0:2]]
    meanx = np.average(xphasepoints)
    varx = np.average(xphasepoints2) - meanx ** 2
    meanerror = np.sqrt(varx / len(xphasepoints))
    print(meanx, meanerror, varx)
    print(max(xphasepoints), min(xphasepoints))
    phaseerr = np.sqrt(meanx ** 2 * w_error ** 2 + w ** 2 * meanerror ** 2)
    print(meanx * w, phaseerr, meanerror * w)
    Dphase = (w * 8.22 ** 2) / (2 * (w * meanx) ** 2)
    Derr = np.sqrt(
        (((8.22) / (w * (meanx ** 2))) ** 2) * (0.05) ** 2 + (((8.22) ** 2) / (2 * (meanx ** 2) * w ** 2)) ** 2 * (
                    (w_error) ** 2) + ((8.22) ** 2 / (w * meanx ** 3)) ** 2 * meanerror ** 2)
    errpercent = (Derr / Dphase) * 100
    print(Dphase, Derr, errpercent)

if T == 360 or T == 480:
    xphasepoints = [x / 10 for x in xval[0:3]]
    print(xphasepoints)
    xphasepoints2 = [(x / 10) ** 2 for x in xval[0:3]]
    meanx = np.average(xphasepoints)
    varx = np.average(xphasepoints2) - meanx ** 2
    meanerror = np.sqrt(varx / len(xphasepoints))
    # print(meanx,meanerror,varx)
    # print(max(xphasepoints), min(xphasepoints))
    phaseerr = np.sqrt(meanx ** 2 * w_error ** 2 + w ** 2 * meanerror ** 2)
    # print(meanx * w, phaseerr, meanerror * w)
    Dphase = (w * (8.22) ** 2) / (2 * (w * meanx) ** 2)
    Derr = np.sqrt(
        (((8.22) / (w * (meanx ** 2))) ** 2) * (0.05) ** 2 + (((8.22) ** 2) / (2 * (meanx ** 2) * w ** 2)) ** 2 * (
                    (w_error) ** 2) + ((8.22) ** 2 / (w * meanx ** 3)) ** 2 * meanerror ** 2)
    errpercent = (Derr / Dphase) * 100
    # print (Dphase,Derr,errpercent)
# standard error of the mean for error in phase and then error in D?
# need to include errors from w, thickness and phase
# print(xval)
yolo = [0.088370, 0.013070517]
# Dphase = [0.12699,0.09052,0.088370,0.0888,0.0927,0.09265]
Dphase = [0.09083, 0.09094, 0.0899, 0.0988, 0.09446, 0.09052, 0.08837, 0.0888, 0.0927, 0.09265, 0.12947, 0.12699,
          0.14778]

# Dphaseerrors = [0.020484943,0.013313071,0.013070517,0.015712981,0.018683137,0.020860171]
Dphaseerrors = [0.00691, 0.00889, 0.00237, 0.002, 0.00186, 0.01354, 0.00668, 0.00348, 0.00257, 0.00206, 0.00304,
                0.00204]

# fractlist=[]
# for x in range(len(Dphase)):
#     fraction = Dphaseerrors[x]/Dphase[x]
#     fractlist.append((1/36)*Dphaseerrors[x]**2)
# #Davgerr = np.sqrt(sum(fractlist))


meanD = np.average(Dphase)
Dphase2 = [x ** 2 for x in Dphase[:]]
varD = np.average(Dphase2) - meanD ** 2
meanDerror = np.sqrt(varD)
meanDerrpercent = meanDerror / meanD * 100
# print(meanD, meanDerror,meanDerrpercent)


# STILL NEED TO COMBINE FOURIER METHOD
alist = []
blist = []
N = 50
for n in range(1, N):
    a_n = spi.simps((2 / T) * np.array(newyv) * np.cos(2 * np.pi * n * np.array(xc) / T), xc)
    b_n = spi.simps((2 / T) * np.array(newyv) * np.sin(2 * np.pi * n * np.array(xc) / T), xc)
    alist.append(a_n)
    blist.append(b_n)

der_wrt_theta_a = spi.simps((2 / T) * np.cos(2 * np.pi * 1 * np.array(xc) / T), xc)
der_wrt_T_a = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.cos(2 * np.pi * 1 * np.array(xc) / T) + (
            4 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.sin(2 * np.pi * 1 * np.array(xc) / T), xc)
der_wrt_T_a3 = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.cos(2 * np.pi * 3 * np.array(xc) / T) + (
            12 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.sin(2 * np.pi * 3 * np.array(xc) / T), xc)
der_wrt_T_a5 = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.cos(2 * np.pi * 5 * np.array(xc) / T) + (
            20 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.sin(2 * np.pi * 5 * np.array(xc) / T), xc)
der_wrt_T_b = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.sin(2 * np.pi * 1 * np.array(xc) / T) - (
            4 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.cos(2 * np.pi * 1 * np.array(xc) / T), xc)
der_wrt_T_b3 = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.sin(2 * np.pi * 3 * np.array(xc) / T) - (
            12 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.cos(2 * np.pi * 3 * np.array(xc) / T), xc)
der_wrt_T_b5 = spi.simps(-(2 / T ** 2) * np.array(newyv) * np.sin(2 * np.pi * 5 * np.array(xc) / T) - (
            20 * np.pi * np.array(xc) / T ** 3) * np.array(newyv) * np.cos(2 * np.pi * 5 * np.array(xc) / T), xc)
a_0_error = np.abs(der_wrt_T_a) * 2
b_0_error = np.abs(der_wrt_T_b) * 2
a_3_error = np.abs(der_wrt_T_a3) * 2
b_3_error = np.abs(der_wrt_T_b3) * 2
a_5_error = np.abs(der_wrt_T_a5) * 2
b_5_error = np.abs(der_wrt_T_b5) * 2
# print(a_0_error,b_0_error)
print(a_3_error, b_3_error)
print(a_5_error, b_5_error)
print(alist)
print(blist)
sum = [np.sum(
    [alist[i] * np.cos((2 * (i + 1) * np.pi * t) / T) + blist[i] * np.sin((2 * (i + 1) * np.pi * t) / T) for i in
     range(N - 1)]) for t in xc]
plt.plot(xc, sum, 'b-')
fourierphase = -np.arctan(alist[0] / blist[0])
fourierphase3 = -np.arctan(alist[2] / blist[2])
fourierphase5 = -np.arctan(alist[3] / blist[3])
if fourierphase < 0:
    fourierphase += np.pi
fourieramp = np.sqrt(alist[0] ** 2 + blist[0] ** 2)
foureramp_err = np.sqrt((alist[0] * (alist[0] ** 2 + blist[0] ** 2) ** -0.5) ** 2 * a_0_error ** 2 + (
            blist[0] * (alist[0] ** 2 + blist[0] ** 2) ** -0.5) ** 2 * b_0_error ** 2)
foureramp_err3 = np.sqrt((alist[2] * (alist[2] ** 2 + blist[2] ** 2) ** -0.5) ** 2 * a_3_error ** 2 + (
            blist[2] * (alist[2] ** 2 + blist[2] ** 2) ** -0.5) ** 2 * b_3_error ** 2)
foureramp_err5 = np.sqrt((alist[4] * (alist[4] ** 2 + blist[4] ** 2) ** -0.5) ** 2 * a_5_error ** 2 + (
            blist[4] * (alist[4] ** 2 + blist[4] ** 2) ** -0.5) ** 2 * b_5_error ** 2)
fourerphase_err = np.sqrt((((1 + (alist[0] / blist[0]) ** 2) ** -0.5) / blist[0]) ** 2 * a_0_error ** 2 + (
            ((1 + (alist[0] / blist[0]) ** 2) ** -0.5) * alist[0] / (blist[0] ** 2)) ** 2 * b_0_error ** 2)
fourerphase_err3 = np.sqrt((((1 + (alist[2] / blist[2]) ** 2) ** -0.5) / blist[2]) ** 2 * a_3_error ** 2 + (
            ((1 + (alist[2] / blist[2]) ** 2) ** -0.5) * alist[2] / (blist[2] ** 2)) ** 2 * b_3_error ** 2)
fourerphase_err5 = np.sqrt((((1 + (alist[4] / blist[4]) ** 2) ** -0.5) / blist[4]) ** 2 * a_5_error ** 2 + (
            ((1 + (alist[4] / blist[4]) ** 2) ** -0.5) * alist[4] / (blist[4] ** 2)) ** 2 * b_5_error ** 2)
print(foureramp_err3, fourerphase_err3)
print(foureramp_err5, fourerphase_err5)
# print(np.sqrt(alist[0]**2+blist[0]**2),fourierphase)

# print(foureramp_err,fourerphase_err)
transfactor1 = fourieramp / 63.7
transfactor3 = np.sqrt(alist[2] ** 2 + blist[2] ** 2) / 21.2
print(transfactor3)
transfactor5 = np.sqrt(alist[4] ** 2 + blist[4] ** 2) / 12.7
# print(damp(transfactor1),dp(fourierphase))
w = 3 * 2 * np.pi / (T)
print(damp(transfactor3), dp(fourierphase3))
w = 5 * 2 * np.pi / (T)
print(damp(transfactor5), dp(fourierphase5))
print(np.sqrt(alist[2] ** 2 + blist[2] ** 2), fourierphase3)
print(np.sqrt(alist[4] ** 2 + blist[4] ** 2), fourierphase5)

# def theta(trans, phase,w,t):
#     return trans * np.sin(phase-w*t)
#
# def anintegrand(n,xc,y):
#     return (2/T)*(yfit*np.cos(((2*np.pi*n*xc)/T)))
#
# def bnintegrand(n,xc,y):
#     return (2/T)*(yfit*np.sin(((2*np.pi*n*xc)/T)))
#
# y1 = anintegrand(1,xc,yfit)
# y2 = bnintegrand(1,xc,yfit)
#
# def an(xc):
#     return spi.simps(y1,xc)
#
# def bn(xc):
#     return spi.simps(y2,xc)
#
# def amp(xc):
#     return np.sqrt(an(xc)**2+bn(xc)**2)
#
# def phase(xc):
#     return -np.arctan(an(xc)/bn(xc))
#
# print (an(xc))
# print(bn(xc))
# print (amp(xc))
# print (phase(xc))


plt.legend(handles=[dataplot, fitplot], loc=1)
plt.show()
