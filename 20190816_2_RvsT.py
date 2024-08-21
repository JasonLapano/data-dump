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
filename = 'RvsT_cool.dat'#raw data from PPMS
thickness = 1#thickness in centimeters
geo = 4.53236 #thickness in centimeters Van der Pauw: pi/ln(2) =4.53236; Hall bar = lengh/Area
#################################
#################################
#################################
#################################
#################################



data = np.genfromtxt(filename, # Specify filename
                     delimiter=',', # How is the data delimited?
                     skip_header=32, # Skip the first 32 rows
                     usecols=(3,4,19,19,20), # Only read these columns
                     comments='Measurement') # Ignore bad lines

data[:,3]=data[:,3]*153/139.5
Rave=geo*(data[:,2]+data[:,3])/2
dRave=np.diff(Rave, n=1)
dRave=np.insert(dRave,0, dRave[0])
ddRave=np.diff(dRave, n=1)
ddRave=np.insert(ddRave,0, ddRave[0])



final_data = {}
final_data=np.stack((data[:,0],data[:,1],data[:,2],data[:,3],Rave,Rave*thickness, dRave, ddRave),axis=1)


plt.plot()
plt.plot(data[:,0],data[:,2])
plt.plot(data[:,0],data[:,3])
plt.plot(data[:,0],Rave)
plt.plot(data[:,0],Rave/geo,'o')
plt.show()


# Export each piece of the split data
head = 'T(K),B(T),Rxx,Ryy,R-sqr,rho,1st-der, 2nd-der'
name = '{:}_averaged.csv'.format(filename)
np.savetxt(name,final_data,delimiter=',', header=head)  
