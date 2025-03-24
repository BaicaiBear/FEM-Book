#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convergence analysis for 2L and 3Q bar element using the bar under compression
shown in Figure 5.13 in Fish's textbook. 

Plot the element length - L2/energy norm cuvers in logarithm scale for both 
linear and quadratic elements, and obtain their convergence rates and the 
intercepts by linear regression.

Created on Thu Apr 30 21:05:47 2020

@author: xzhang
"""

from Bar1D import FERun
from Exact import ErrorNorm_HeatConductionBar
import FEData as model
from Bar1DElem import Nmatrix1D, Bmatrix1D
import numpy as np
import matplotlib.pyplot as plt

# Json data files for 2L element
files_2L = ("4-elements.json", "8-elements.json",
            "16-elements.json")    


# Run FE analysis for all files using 2L element
n2L = len(files_2L)
h2 = np.zeros(n2L)
L2Norm2 = np.zeros(n2L)
FluxDifference2 = np.zeros(n2L)
for i in range(n2L):
    FERun("Convergence/HeatConductionBar/"+files_2L[i])

    # Calculate error norms for convergence study
    h2[i], L2Norm2[i] = ErrorNorm_HeatConductionBar()

    # Calculate the difference of heat flux at the right end of the bar
    # compute the heat flux at the bar's right end in a manner similar to stress calculation
    # here we take the last element, evaluate flux at x = xe[-1], and compare with exact
    e = model.nel - 1
    de = model.d[model.LM[:,e] - 1]
    IENe = model.IEN[:,e] - 1
    xe = model.x[IENe]

    N  = Nmatrix1D(xe[-1], xe)
    B  = Bmatrix1D(xe[-1], xe)
    kE = N @ model.E[IENe]      # thermal conductivity at the element's right end

    flux_FE = -kE * (B @ de)    # FE heat flux (minus sign for q = -k dT/dx)
    flux_exact = 0  # user-defined exact solution

    FluxDifference2[i] = flux_FE - flux_exact


log_h2 = np.log(h2)
log_err2 = np.log(L2Norm2)

plt.figure()
plt.plot(log_h2, log_err2, 'o-', label='Linear Element')
plt.xlabel('log(h)')
plt.ylabel('log(L2 Error)')
plt.title('Linear Element Error Convergence')
slope2, _ = np.polyfit(log_h2, log_err2, 1)
plt.plot(log_h2, slope2*log_h2, '--', label=f'Slope={slope2:.2f}')
plt.legend()
plt.grid(True)
plt.show()

# Plot flux difference in log-log scale (absolute value in case of negative flux)
plt.figure()
plt.plot(h2, FluxDifference2, 's-', label='Boundary Error')
plt.xlabel('h')
plt.ylabel('Flux Error')
plt.title('Boundary Error Convergence')
slope_flux, _ = np.polyfit(h2, FluxDifference2, 1)
plt.plot(h2, slope_flux*h2, '--', label=f'Flux Slope={slope_flux:.2f}')
plt.legend()
plt.grid(True)
plt.show()

# Linear regression 
print("The L2/energy error norms are ")


a, C = np.polyfit(np.log(h2),np.log(L2Norm2),1)
print("    Linear element    : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(FluxDifference2),1)
print("    Linear element    : e_n = %e h^%g" %(np.e**C, a))


# Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
# which can be added into your LaTex source code by "\input{fe_plot.tex}"
import tikzplotlib
tikzplotlib.save("fe_convergence.tex")


# Print error norms obtained by the linear element and quadratic element
#    with different element size
print("\nError norms of linear elements")
print('%13s %13s %13s' %('h','L2Norm','FluxDifference'))
for i in range(len(h2)):
    print('%13.6E %13.6E %13.6E' %(h2[i], L2Norm2[i], FluxDifference2[i]))
