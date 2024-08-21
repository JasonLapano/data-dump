# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 11:06:21 2021

@author: 5ul
"""


import function_list
import QuantumOscillation_fitting
import importer
import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy 
import scipy as sp

#from scipy.misc import derivative
#plt.style.use('custom')


#head = 'T(K),B(T),1/B, Rxy(ohms),Rxy-background(ohm), rho(ohm*cm), Magneto,dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw'
# S:
#     0 = T
#     1 = B
#     2 = 1/B
#     3 = Rxy
#     4 = Rxy-AH
#     5 = rho
#     6 = MR
#     7 = dRxy/dB
#     8 = dRsquare/dB
#     9 = drho/dB
#     10 = Rave
#     11 = Rsym
#     12 = Rxy_raw
#     13 = Rxx_raw
#     14 = Ryy_raw
#     15 = delta_Rxx
#     16 = delta_Rxy
#     17 = Rxx_ave
#     18 = Rxy_ave



#C:\Users\5ul\Google Drive\RuCl3_graphite\RuCl3_Graphite\20210803_1_stg2\dil fridge\drive-download-20210809T180252Z-001\00000_500mK.da
#funf.RvsH_PM_import(filename, cols=[5,6,3, 3, 4], head=2, geo=1, deli='\t', Bmax=14.01)

filename = r'D:\PPMS\RuCl3_Graphite\20210917\RvsH_2K.dat'
s = importer.RvsH_FM_import(filename, cols=[0,1,3,3,3], head=1, geo=1, deli='\t', Bmax=9.1)
#%%

for i in range(len(s)):
    sBin = function_list.bin_fun(s[i], 1)
    x = np.size(sBin) #gets size of s for trimming 
    
    #Fits the parabolic background
    params, params_covariance, fit, iB, subs = QuantumOscillation_fitting.SDH_osc(sBin[1:(x-1), :])
    deg = str(np.round(s[i][0,0]))
    
    
    #interpolates the curve with respect to 1/B
    minfield = 0.3 #where oscillations start
    maxfield = 9.1    # where oscillations stop or max field value
    res = 100000 # resolution of fitting
    x = np.linspace(1/maxfield, 1/minfield, num=res, endpoint=True)
    f3 = scipy.interpolate.interp1d(iB[abs(iB)<=5], subs[abs(iB)<=5], kind='linear')
    # plt.figure()
    # plt.plot(x, f3(x))
    
    l = np.size(iB)
    ex = np.zeros((res,2))
    ex[:, 0] = x
    ex[:, 1] = f3(x)
    
    # Put the temperature in the filename
    name = filename.strip('.dat')+'_{:}_cond.csv'.format(np.round(s[i][1,0]))
    #Save the file
    np.savetxt(name, ex, delimiter=',')  
    
    plt.figure(1)
    plt.plot(x, f3(x)+i*.00000001, label=deg)
    
    
    plt.figure(3)
    fy = sp.fft(f3(x))
    max_freq = (1/minfield- 1/maxfield)/res
    fx = np.fft.fftfreq(res, (1/minfield-1/maxfield)/res)
    plt.semilogy(fx, abs(fy)+i*.0001, label=deg)
    plt.xlim(0, 10000) 
  
#%%
    
plt.rcParams.update({'font.size': 12})   
for i in range(len(s)):

    plt.figure(1)
    plt.plot(x, f3(x)+i*.001, label=deg)
    
    plt.figure(2)
    fy = sp.fft(f3(x))
    plt.plot(abs(fy)+i*.001, label=deg)

plt.figure(2)    
plt.rcParams.update({'font.size': 12})   
plt.xlim(0.11, 0.2)   

#%% 
print(np.size(s[0][0,:]))

#%%
for i in range(len(s)):
    
    plt.figure(1)
    plt.plot(s[i][:,1],s[i][:,17]+i*.0000002)
    plt.xlim([-9,9])
    
    plt.figure(1)
    


