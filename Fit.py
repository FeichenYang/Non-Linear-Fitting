#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

dataset = pd.read_csv("Creep for cartilage.csv")

time = dataset.iloc[:,0].values
strain_true = dataset.iloc[:,1].values

"""
the zero point of the time is defined as the time when the applied pressure 
reaches the applied stress 0.04 MPa
the baseline thickness was defined as the thickness at near zero pressure (the 
samples were preloaded with 0.02N force for 15 minutes before tests)
"""

#sample height 3 mm
h = 0.00285

from scipy.optimize import leastsq

#calculate the function, summing the first 5 items
#input argument: p is the parameters, or the constants; x is the variables.
def strain_fit (p,x):
    #H_A is aggregate modulus, tau is characteristic time
    H_A, tau = p
    time = x
    u_over_h = 0.04/H_A
    #the equation is u/h times that really long thing
    factor = 1
    #sum over n = 0 to 5
    for i in range(11):
        factor -= 2*(sp.pi**-2*(i+0.5)**-2*sp.exp(-(sp.pi**2*(i+0.5)**2)*time/tau))
    return u_over_h*factor

#calculate error, input argument: 
def error(p,x,y):
    return strain_fit(p,x)-y

# initial parameter with Hattie's results: 0.71 MPa H_A, 11.6 min characteristic time
p0 = [0.71, 696]

#non-linear fit, least square
p_fit = leastsq(error, p0, args = (time, strain_true))

#plot the fitted line and compare with original line
strain_fitted = strain_fit(p_fit[0], time)

plt.plot(time, strain_fitted,label = "fitted line", color = "red")
plt.plot(time, strain_true, label = "original line", color = "black")
plt.legend()
plt.show()

# output measurements
tau = p_fit[0][1]
Aggregate_Modulus = p_fit[0][0]

permeability = h**2/(Aggregate_Modulus*1E6)/tau
