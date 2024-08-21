# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 09:11:42 2019

@author: m4b
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
filename = '20210302_RvsT_325-2_0_T.dat'#raw data from PPMS
direction = -1 #-1 takes -B->+B, +1 takes +B->-B
Bmax = 9.001 #max field used; important for interpolation fucntion
thickness = 1*10**-7#thickness in centimeters
geo = 1 #thickness in centimeters Van der Pauw: pi/ln(2) =4.53236; Hall bar = lengh/Area
saturation = 20 # % of data used to fit the linear background (i.e. for a 10T scan, 10 would give a fit over 9-10T data range)




scans = 0 #  Set to 1 if you are running a series of temperature scans
config = 1 #set to 1, 2 or 3 to choose what channels are being read for the data (see below)

#################################
#################################
#################################
#################################
#################################

if config == 1:
    
    ## assuming with lock-in, there is only Rxx and Rxy. Time is in the last column, may have to work around later
    data = np.genfromtxt(filename, # Specify filename
                     delimiter='\t', # How is the data delimited?
                     skip_header=1, # Skip the first 32 rows
                     usecols=(1,2,8,8,9,0), # Only read these columns
                     comments='Measurement') # Ignore bad lines

# if config == 2:
#     data = np.genfromtxt(filename, # Specify filename
#                      delimiter=',', # How is the data delimited?
#                      skip_header=32, # Skip the first 32 rows
#                      usecols=(3,4,21,21,22), # Only read these columns
#                      comments='Measurement') # Ignore bad lines
    
# if config == 3:
#     data = np.genfromtxt(filename, # Specify filename
#                      delimiter=',', # How is the data delimited?
#                      skip_header=32, # Skip the first 32 rows
#                      usecols=(3,4,19,20,21), # Only read these columns
#                      comments='Measurement') # Ignore bad lines

data[:,3]=data[:,3]
Rave=geo*(data[:,2]+data[:,3])/2
dRave=np.diff(Rave, n=1)
dRave=np.insert(dRave,0, dRave[0])
ddRave=np.diff(dRave, n=1)
ddRave=np.insert(ddRave,0, ddRave[0])



final_data = {}
final_data=np.stack((data[:,0],data[:,1],data[:,2],data[:,3],Rave,Rave*thickness, dRave, ddRave),axis=1)

plt.rcParams.update({'font.size': 10})
ax1 = plt.subplot(211)
plt.plot(data[:,0],data[:,2])
plt.plot(data[:,0],data[:,3])
plt.plot(data[:,0],Rave/geo)
ax1.set_title("Raw Data")
ax1.set_xlabel("Temperature (K))")
ax1.set_ylabel("Resistance(Ohm)")
ax2 = plt.subplot(212)
plt.plot(data[:,0],Rave)
ax2.set_title("Symmetrized w/ Geo Factor")
ax2.set_xlabel("Temperature (K))")
ax2.set_ylabel("Resistance(Ohm)")
plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
plt.show()
plt.savefig('RvsT_scan.png', dpi = 400)


# Export each piece of the split data
head = 'T(K),B(T),Rxx,Ryy,R-sqr,rho,1st-der, 2nd-der'
name = '{:}_averaged.csv'.format(filename)
np.savetxt(name,final_data,delimiter=',', header=head)  
