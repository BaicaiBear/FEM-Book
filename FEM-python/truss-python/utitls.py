#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides utilities used by FE analysis.
  1. assembly: Global stiffness matrix assembly.
  2. solvedr: Solving the stiffness equations by the reduction approach.

Created on Sat May 9 17:39:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import FEData as model


def assembly(e, ke):
    """
    Assemble element stiffness matrix.
    
    Args:
        e   : (int) Element number
        ke  : (numpy(nen*ndof,nen*ndof)) element stiffness matrix
    """
    model.K[np.ix_(model.LM[:,e], model.LM[:,e])] += ke

def solvedr():
    """
    Partition and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd = model.nd; neq=model.neq
    K_E = model.K[0:nd, 0:nd]
    K_F = model.K[nd:neq, nd:neq]
    K_EF =model. K[0:nd, nd:neq]
    f_F = model.f[nd:neq]
    d_E = model.d[0:nd]
    
    # solve for d_F
    d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E) 

    # reconstruct the global displacement d
    model.d = np.append(d_E,d_F)
    
    # compute the reaction r
    f_E = K_E@d_E + K_EF@d_F
    
    # write to the workspace
    print('\nsolution d');  print(model.d)
    print('\nreaction f =', f_E)
    
    return f_E

def solvepf(beta=1e16):
    """
    Set punishment function and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd = model.nd; neq=model.neq
    K = model.K.copy()
    f = model.f.copy()
    for i in range(nd):
        K[i,i] = beta
        f[i] = beta*model.d[i]
    # a node cannot have force and displacement confinement at the same time
    # solve for d
    d = np.linalg.solve(K, f)
    model.d = d
    # compute the reaction r
    f = K@d
    # write to the workspace
    print('\nsolution d');  print(model.d.T)
    print('\nreaction f =', f)

    return f

def solvepf_3D(beta=1e16):
    """
    Set punishment function and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd = model.nd; neq=model.neq
    K = model.K.copy()
    f = model.f.copy()
    for i in range(nd):
        K[i,i] = beta
        f[i] = beta*model.d[i]
    for j in range(neq):
        if(j%3==2):
            K[j,j] = beta
            f[j] = 0
    # a node cannot have force and displacement confinement at the same time
    # solve for d
    d = np.linalg.solve(K, f)
    model.d = d
    # compute the reaction r
    f = K@d
    # write to the workspace
    print('\nsolution d');  print(model.d.T)
    print('\nreaction f =', f)

    return f