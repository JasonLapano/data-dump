# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter

plt.close('all')
#put all filenames of scans here.  Filepath in 1st column, name in second as .xy file, 3rd as name for plotting.  filename should start with \\
filenames = np.array([[r'C:\Users\5ul\Google Drive\MnTe\XRD\20210504_2', r'\20210504_2.xy', '0504_2-S'],
                      [r'C:\Users\5ul\Google Drive\MnTe\XRD\20210505_1', r'\20210505_1_u00_1.xy', '0505_1-S'], 
                      [r'C:\Users\5ul\Google Drive\MnTe\XRD\20210506_1', r'\20210506_1_u00.xy', '0506_1-S'], 
                      [r'C:\Users\5ul\Google Drive\MnTe\XRD\20210429_1', r'\5to105std.xy', '0429_1-Fe-S'],
                      [r'C:\Users\5ul\Google Drive\MnTe\XRD\20210429_2', r'\5to105std.xy', '0429-2-Fe-MPMS'], 
                      [r'C:\Users\5ul\Google Drive\MnTe\XRD\20210514_2', r'\202105014_2_u00.xy', '0514-1-Fe']])

fig1, ax1 = plt.subplots()
i = 0;
for filename in filenames:
    
    name = filename[0] + filename[1]
    twotheta, intensity = np.loadtxt(name, delimiter=' ', unpack=True)
    intensity = (intensity + 2) *10**-i
    i = i + 2
    plt.plot(twotheta, intensity, label = 'XRD Scan', linewidth=0.5)
    plt.yscale('log')
    plt.xlabel('2Theta')
    plt.ylabel('Intensity')
    plt.ion()
   

plt.legend(filenames[:, 2], prop={'size':10})
ax1.xaxis.set_major_locator(MultipleLocator(5))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(1))  
ax1.tick_params(direction='in', which='both', right=True, top=True, labelright=False, labeltop=False)
ax1.plot(label = '2theta' )

ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, numticks=10))  
#ax1.tick_params(direction='in', which='both', right=True, top=True, labelright=False, labeltop=False)
plt.xlim(20, 60)
plt.show()



a = 12.994
l = 2 * a * np.sin(twotheta*np.pi/180/2) / (1.5401)
fnl = np.zeros((l.size, 5))

fnl[:, 2] = l
fnl[:, 3] = intensity

plt.savefig(r'C:\Users\5ul\Google Drive\MnTe\scans.png')
plt.savefig(r'C:\Users\5ul\Google Drive\MnTe\scans.svg')
