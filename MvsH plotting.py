# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:14:08 2019

@author: Jason
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 09:16:52 2017

@author: Matthew Brahlek_updated 9/26/2016
"""
import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
#from scipy.misc import derivative
#plt.style.use('custom')

plt.close("all")


#################################
#################################
#######input information#########
#################################
#################################
filenames =['SbMBT_25p@5K.dat']#raw data from PPMS

direction = -1 #-1 takes +B->-B, +1 takes -B->+B
Bmax = 7.01 #max field used; important for interpolation fucntion
thickness =  10 #thickness in nanometers
area = 0.25 #sample surface area in centimeters
saturation = 50 # % of data used to subtract the linear background
#################################
#################################
#################################
#################################
#################################


volume = abs(np.cos(60))*.436**2*2.36*10**-21 #unit cell volume in cm
atoms = 1 #atoms/unit cell
uB = 1.07*10**20 #convert EMU to uB
#MvsH as as function of T

#MvsH as as function of T
fig, ax1 = plt.subplots()
plt.title('MvsH')
plt.xlabel('Field (Tesla)')
plt.ylabel('Magnetization')
fig, ax2 = plt.subplots()
plt.title('MvsH substrate substracted')
plt.xlabel('Field (Tesla)')
plt.ylabel('Magnetization')
for filename in filenames:
    data = np.genfromtxt(filename, # Specify filename
                     delimiter=',', # How is the data delimited?
                     skip_header=29, # Skip the first 32 rows
                     usecols=(2,3,38,39), # Only read these columns
                     comments='"') # Ignore bad lines
    
    plt.figure(1)   
    plt.rcParams.update({'font.size': 15})
    data[:, 2] = data[:, 2]/(area*thickness*10**-7); #EMU/cm^3
    data[:, 2] = data[:, 2]*volume/atoms*uB
    #plot raw MvsH
    temp = np.array2string(np.round_(data[0, 0], decimals=1))
    ax1.plot(data[:, 1]*10**-4, data[:, 2], label = (temp+' K'))
    
    plt.figure(2)    
    linear = np.where(data[:, 1]>np.amax(data[:, 1])*(1-saturation*.01)) #finds where data is saturated for linear fit
    linear = linear[0] #converts from tuple to float
    #substract linear background and plot
    z = np.polyfit(data[linear, 1], data[linear, 2], 1);
    subs = data[:, 2] - data[:, 1]*z[0] 
    ax2.plot(data[:, 1]*10**-4, subs, label = (temp+' K'))
    
    # creates name
    dataname = filename.strip('.dat')+'.csv'
    rdimage = filename.strip('.dat')+'_RD'+'.png'
    bsimage = filename.strip('.dat')+'_BS'+'.png'
    
    
    leg1 = ax1.legend(); 
    leg2 = ax2.legend();   
    
    plt.figure(1)
    plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
    plt.savefig(rdimage, dpi = 400)
    
    plt.figure(2)
    plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
    plt.savefig(bsimage, dpi = 400)    
    
    export=np.empty((len(data), 4 ))
    export[:,0] = data[:, 1]*10**-4
    export[:,1] = data[:, 2]
    export[:,2] = subs
    export[:,3] = data[:, 1]*z[0]
    
    head = 'Field(B), Magnetization (uB), Substracted (uB), Background (uB)'
    np.savetxt(dataname, export, delimiter=',', header=head)
   
    