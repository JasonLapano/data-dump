# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.close('all')
#put all filenames of scans here.  Filepath in 1st column, name in second as .xy file, 3rd as name for plotting.  filename should start with \\
filenames = np.array([[r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200313_2_20200314_1\20200314_1_10%', r'\RvsH_2K_output_2.5K.csv', '10%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200420_1_12p', r'\20190420_1_RvsH_2.5K_output_2.5K.csv', '12%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200325_2_15p', r'\RvsH_2k_edi_output_2.5K.csv', '15%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200313_2_20200314_1\30300313_2_20%', r'\RvsH_2K_output_2.5K.csv', '20%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200326_1_25p', r'\RvsH_2k_output_2.0K.csv', '25%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200314_2_20200315_1\20200314_2', r'\RvsH_2kb_output_2.5K.csv', '30%'],
                      [r'C:\Users\5ul\Google Drive\Sb-BMT\Transport\20200314_2_20200315_1\20200315_1', r'\RvsH_2kb_output_2.5K.csv', '40%']])

#fig1= plt.subplots(nrows=3, ncols=filenames.shape[0])
matplotlib.rcParams.update({'font.size': 8})
i = 0;
for filename in filenames:
    
    name = filename[0] + filename[1]
    Temp, Field, Rxy, rho, dRxy, dRsqr, drho, Rave, Rsym, Rxy_raw, Rxx_raw, Ryy_raw = np.loadtxt(name, delimiter=',', unpack=True)
    ax1=plt.subplot(3, filenames.shape[0], i+1)
    plt.plot(Field, Rsym)
    ax1.set_title(filename[2])
    ax1.set_xlabel("Field(B)")
    ax1.set_ylabel("Resistance(Ohm)")
    ax1=plt.subplot(3, filenames.shape[0], filenames.shape[0]+i+1)
    plt.plot(Field, Rxy)
#    ax1.set_title("Hall Resistance")
#    ax1.set_xlabel("Field(B)")
#    ax1.set_ylabel("Resistance(Ohm)")
#    ax1.set_ylim([-2000, 2000])
    ax1=plt.subplot(3, filenames.shape[0], 2*filenames.shape[0]+i+1)
    plt.plot(Field, Rxy)
#    ax1.set_title("Lowfield Hall")
#    ax1.set_xlabel("Field(B)")
#    ax1.set_ylabel("Resistance(Ohm)")
    ax1.set_xlim([-3, 3])
    ax1.set_ylim([-1000, 1000])
    
    i = i+1
#    plt.plot(twotheta, intensity, label = 'XRD Scan', linewidth=0.5)
#    plt.yscale('log')
#    plt.xlabel('2Theta')
#    plt.ylabel('Intensity')
#    plt.ion()
plt.tight_layout(pad=.2, h_pad=2, w_pad=None, rect=None)

#   
#
#plt.legend(filenames[:, 2], prop={'size':10})
#ax1.xaxis.set_major_locator(MultipleLocator(5))
#ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
#ax1.xaxis.set_minor_locator(MultipleLocator(1))  
#ax1.tick_params(direction='in', which='both', right=True, top=True, labelright=False, labeltop=False)
#ax1.plot(label = '2theta' )
#
#ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
##ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
#ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, numticks=10))  
##ax1.tick_params(direction='in', which='both', right=True, top=True, labelright=False, labeltop=False)
#plt.xlim(0, 80)
#plt.show()
#
#
#
#a = 12.994
#l = 2 * a * np.sin(twotheta*np.pi/180/2) / (1.5401)
#fnl = np.zeros((l.size, 5))
#
#fnl[:, 2] = l
#fnl[:, 3] = intensity
#
plt.savefig(r'C:\Users\5ul\Google Drive\Sb-BMT\Figures\Magnetotransport\scans.png')
plt.savefig(r'C:\Users\5ul\Google Drive\Sb-BMT\Figures\Magnetotransport\scans.svg')
