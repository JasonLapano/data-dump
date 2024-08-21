# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:59:07 2019

@author: 5ul
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

#[filelocation, filename, shortname, T]

filenames = [[r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_2.0K.csv', '2K', 2],
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_10.0K.csv', '10K', 10],
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_20.0K.csv', '20K', 20], 
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_50.0K.csv', '50K', 50],
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_100.0K.csv', '100K', 100],
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_200.0K.csv', '200K', 200],
             [r'C:\Users\5ul\Google Drive\Bi2Se3-CoSe-Bi2Se3\transport\20191209_1', r'\RvsH_2K_output_300.0K.csv', '300K', 300]]

for n in range(len(filenames)):
    filename = filenames[i][0]+filenames[i][2]