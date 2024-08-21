"""
# -*- coding: utf-8 -*-
"""
"""
Created on Thu Dec  5 13:24:12 2019

@author: 5ul
"""
"""

# -*- coding: utf-8 -*-
"
Created on Tue Sep 26 09:16:52 2017

@author: Matthew Brahlek_updated 9/26/2016
"""
import csv
import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
#from scipy.misc import derivative
#plt.style.use('custom')
import sympy as sympy
from sympy import symbols

#Define symbols
q, n1, m1, m1, B, n2, m2 = symbols('q n1 m1 m1 b n2 m2' )

#Define conductance
Gxx = q*((-n1*m1)/(1+m1**2*B**2) + (n2+m2)/(1+m2**2*B**2))
Gxy = 1*B*(n1*m1**2/(1-m1**2*B**2) + n2*m2/(1+m2**2*B**2))

#Define resistance in terms of conductance
Rxx = Gxx/(Gxx**2 + Gxy**2)
Rxy = Gxy/(Gxx**2+Gxy**2)

Rxx0 = Rxx.subs(B, 0)
dRxy = sympy.diff(Rxy, B)