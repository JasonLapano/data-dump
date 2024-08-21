# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 13:12:03 2021

@author: 5ul
"""

import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
from datetime import datetime
import time

#lockin data
lockin = pd.read_csv(r'D:\PPMS\RuCl3_Graphite\20210917\jason_real.dat')

#ppms file
ppms = np.genfromtxt(r'D:\PPMS\RuCl3_Graphite\20210917\full_dataset.dat',
                     delimiter=',', # How is the data delimited?
                     skip_header=32, # Skip the first 32 rows
                     usecols=(1,3,4,5,19,19,20), # Only read these columns
                     comments='Measurement')

#%%
# converting time of scan from 0 to the end of the lockin data in seconds
day = []
#finds the first datapoint time
d = datetime.strptime(lockin.iloc[0,0] + lockin.iloc[0,1], '%Y/%m/%d %H:%M:%S.%f')
#ztime = zero time
ztime = time.mktime(d.timetuple())

#this determines the time for the lockin "day".  Must be done row by row
for i in range(len(lockin.iloc[:,3])): # of datapoints
    d = datetime.strptime(lockin.iloc[i,0] + lockin.iloc[i,1], '%Y/%m/%d %H:%M:%S.%f')
    ltime = time.mktime(d.timetuple())-ztime 
    day = np.append(day, ltime) #day = lockintime - ztime
    
    
    
    
#converts ppms data from 0 to the end of the scan in seconds


#%%
day2 = abs(day - min(day) - 21.5) #this takes care of the lag between the data. Play around with the time to add or substract to get the lockins and ppms at the same start point
#negative times throw an error, so use absolute value to just wrap the times back aroudn to positive
ppms[:, 0] = ppms[:,0] - ppms[0,0]
interptime = day # this can be deleted i think
#put in ppms data here

#creates the interpolation function for the values off the ppms
ftemp = interp1d(ppms[:,0], ppms[:,1], kind='cubic')
ffield = interp1d(ppms[:,0], ppms[:,2], kind='cubic')
frotation = interp1d(ppms[:,0], ppms[:,3], kind='cubic')

#creates the interpolated temperature, field and rotation data to match the lockin data density
#the range is fixed as both day2 and ppms times should start at 0, and only interpolates where day2 is < ppms[:,0]
lockin_temp = ftemp(day2[day2<max(ppms[:,0])])
lockin_field = ffield(day2[day2<max(ppms[:,0])])
lockin_rotation = frotation(day2[day2<max(ppms[:,0])])


plt.plot(lockin_field, lockin.iloc[:,4].to_numpy()[day2<max(ppms[:,0])])

#%%
#exports data values
x = lockin_field
y = lockin.iloc[:,2:5].to_numpy()[day2<max(ppms[:,0])]


data = np.zeros((len(lockin_field), 7))
data[:,0] = lockin_temp
data[:,1] = lockin_field
data[:, 2] = lockin_rotation
data[:, 3:6] = lockin.iloc[:,2:5].to_numpy()[day2<max(ppms[:,0])]

head = 'T,Field, Position,X,Y,R,theta'
np.savetxt(r'D:\PPMS\RuCl3_Graphite\20210917\converted_data.txt', data, delimiter=',', header=head)       

#%%


