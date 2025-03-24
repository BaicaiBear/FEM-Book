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
from Exact import ErrorNorm_CompressionBar

import numpy as np
import matplotlib.pyplot as plt

# Json data files for 2L element
files_2L = ("2-elements.json", "4-elements.json", "8-elements.json",
            "16-elements.json", "32-elements.json")    

# Json data files for 3Q element
files_3Q = ("2-elements-3Q.json", "4-elements-3Q.json", "8-elements-3Q.json",
            "16-elements-3Q.json")    

# Run FE analysis for all files using 2L element
n2L = len(files_2L)
h2 = np.zeros(n2L)
L2Norm2 = np.zeros(n2L)
EnNorm2 = np.zeros(n2L)
for i in range(n2L):
    FERun("Convergence/CompressionBar/"+files_2L[i])

    # Calculate error norms for convergence study
    h2[i], L2Norm2[i], EnNorm2[i] = ErrorNorm_CompressionBar()

# Run FE analysis for all files using 3Q element
n3Q = len(files_3Q)
h3 = np.zeros(n3Q)
L2Norm3 = np.zeros(n3Q)
EnNorm3 = np.zeros(n3Q)
for i in range(n3Q):
    FERun("Convergence/CompressionBar/"+files_3Q[i])

    # Calculate error norms for convergence study of the bar under compression
    h3[i], L2Norm3[i], EnNorm3[i] = ErrorNorm_CompressionBar()

fig, axs = plt.subplots(2, 2, figsize=(8, 6), tight_layout=True)

# Linear L2 error
log_h2 = np.log(h2)
log_l2_2 = np.log(L2Norm2)
axs[0,0].plot(log_h2, log_l2_2, 'o-', label='Linear L2')
slope2, _ = np.polyfit(log_h2, log_l2_2, 1)
axs[0,0].plot(log_h2, slope2*log_h2, '--', label=f'Slope={slope2:.2f}')
axs[0,0].set_xlabel('log(h)')
axs[0,0].set_ylabel('log(L2 Error)')
axs[0,0].set_title('Linear L2 Convergence')
axs[0,0].legend()
axs[0,0].grid(True)

# Quadratic L2 error
log_h3 = np.log(h3)
log_l2_3 = np.log(L2Norm3)
axs[0,1].plot(log_h3, log_l2_3, 'o-', label='Quadratic L2')
slope3, _ = np.polyfit(log_h3, log_l2_3, 1)
axs[0,1].plot(log_h3, slope3*log_h3, '--', label=f'Slope={slope3:.2f}')
axs[0,1].set_xlabel('log(h)')
axs[0,1].set_ylabel('log(L2 Error)')
axs[0,1].set_title('Quadratic L2 Convergence')
axs[0,1].legend()
axs[0,1].grid(True)

# Linear energy error
log_en2 = np.log(EnNorm2)
axs[1,0].plot(log_h2, log_en2, 'o-', label='Linear Energy')
slope2e, _ = np.polyfit(log_h2, log_en2, 1)
axs[1,0].plot(log_h2, slope2e*log_h2, '--', label=f'Slope={slope2e:.2f}')
axs[1,0].set_xlabel('log(h)')
axs[1,0].set_ylabel('log(Energy Error)')
axs[1,0].set_title('Linear Energy Convergence')
axs[1,0].legend()
axs[1,0].grid(True)

# Quadratic energy error
log_en3 = np.log(EnNorm3)
axs[1,1].plot(log_h3, log_en3, 'o-', label='Quadratic Energy')
slope3e, _ = np.polyfit(log_h3, log_en3, 1)
axs[1,1].plot(log_h3, slope3e*log_h3, '--', label=f'Slope={slope3e:.2f}')
axs[1,1].set_xlabel('log(h)')
axs[1,1].set_ylabel('log(Energy Error)')
axs[1,1].set_title('Quadratic Energy Convergence')
axs[1,1].legend()
axs[1,1].grid(True)

plt.show()

# Linear regression 
print("The L2/energy error norms are ")


a, C = np.polyfit(np.log(h2),np.log(L2Norm2),1)
print("    Linear element    : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(EnNorm2),1)
print("    Linear element    : ||e||_en = %e h^%g" %(np.e**C, a))

a, C = np.polyfit(np.log(h3),np.log(L2Norm3),1)
print("    Quadratic element : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h3),np.log(EnNorm3),1)
print("    Quadratic element : ||e||_en = %e h^%g\n" %(np.e**C, a))

# Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
# which can be added into your LaTex source code by "\input{fe_plot.tex}"
import tikzplotlib
tikzplotlib.save("fe_convergence.tex")


# Print error norms obtained by the linear element and quadratic element
#    with different element size
print("\nError norms of linear elements")
print('%13s %13s %13s' %('h','L2Norm','EnNorm'))
for i in range(len(h2)):
    print('%13.6E %13.6E %13.6E' %(h2[i], L2Norm2[i], EnNorm2[i]))

print("\nError norms of quadratic elements")
print('%13s %13s %13s' %('h','L2Norm','EnNorm'))
for i in range(len(h3)):
    print('%13.6E %13.6E %13.6E' %(h3[i], L2Norm3[i], EnNorm3[i]))
