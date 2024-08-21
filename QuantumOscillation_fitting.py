# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 15:44:11 2021
Fit Linear Transport data
@author: 5ul
""" 
import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import optimize
from scipy.signal import savgol_filter

def test_func(x, a, b, c, Bcrit, m):
    return a*(x**2) + b*(x**2)**(1/2)*  (((x**2)**0.5)/Bcrit)**m  / (1+ (((x**2)**0.5)/Bcrit)**m)  + c

def SDH_osc(s):
    k = 0  #select which temperature data this is

    HF =  max(s[:, 1])          #High field
    MRsign = 1                  #1 for positive MR, -1 for negative
    
    # linear = np.where(s[:, 1]>np.amax(s[:, 1])*.8) #finds where data is saturated for linear fit
    # linear = linear[0]
    # z = np.polyfit(s[linear, 1], s[linear, 5], 1);
    # subs = s[:, 5] - abs(s[:, 1]*z[0])*MRsign 
    
    #smoothed = savgol_filter(s[:, 15], 11, 3)
    smoothed = s[:, 5]
    params, params_covariance = optimize.curve_fit(test_func, s[:,1], smoothed, maxfev=100000)
    fit = test_func(s[:,1], params[0], params[1], params[2], params[3], params[4])
    subs = smoothed - fit
   # plt.figure()
   # plt.plot(s[:,1], fit)
    #plt.plot(s[:, 1], smoothed)
    
   # plt.figure()
    iB = 1/s[:, 1]
    # plt.plot(iB, subs)
    return params, params_covariance, fit, iB, subs
    