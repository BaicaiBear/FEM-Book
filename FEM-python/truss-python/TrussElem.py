#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to calculate element stiffness matrix

Created on Sat May 9 15:40:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import FEData as model
import numpy as np

def TrussElem(e):
    '''
    calculate element stiffness matrix

    Args:
        e : (int) element number

    Returns:
        ke : (numpy(nen,nen)) element stiffness matrix
    '''
    # constant coefficient for each truss element
    const = model.CArea[e]*model.E[e]/model.leng[e]

    # calculate element stiffness matrix
    if model.ndof == 1:
        ke = const*np.array([[1, -1], [-1, 1]])
    elif model.ndof == 2:
        IENe = model.IEN[e] - 1
        xe = model.x[IENe]
        ye = model.y[IENe]
        s = (ye[1] - ye[0])/model.leng[e]
        c = (xe[1] - xe[0])/model.leng[e]

        s_s = s*s
        c_c = c*c
        c_s = c*s

        ke = const*np.array([[c_c, c_s, -c_c, -c_s],
                             [c_s, s_s, -c_s, -s_s],
                             [-c_c, -c_s, c_c, c_s],
                             [-c_s, -s_s, c_s, s_s]])
    elif model.ndof == 3:
        IENe = model.IEN[e] - 1
        xe = model.x[IENe]
        ye = model.y[IENe]
        ze = model.z[IENe]
        l = model.leng[e]
        s_x = (xe[1] - xe[0]) / l
        s_y = (ye[1] - ye[0]) / l
        s_z = (ze[1] - ze[0]) / l

        s_xx = s_x * s_x
        s_yy = s_y * s_y
        s_zz = s_z * s_z
        s_xy = s_x * s_y
        s_xz = s_x * s_z
        s_yz = s_y * s_z

        ke = const * np.array([[s_xx, s_xy, s_xz, -s_xx, -s_xy, -s_xz],
                               [s_xy, s_yy, s_yz, -s_xy, -s_yy, -s_yz],
                               [s_xz, s_yz, s_zz, -s_xz, -s_yz, -s_zz],
                               [-s_xx, -s_xy, -s_xz, s_xx, s_xy, s_xz],
                               [-s_xy, -s_yy, -s_yz, s_xy, s_yy, s_yz],
                               [-s_xz, -s_yz, -s_zz, s_xz, s_yz, s_zz]])
        # insert your code here for 3D
        # ...
        # pass # delete or comment this line after your implementation for 3D
    else:
        raise ValueError("The dimension (ndof = {0}) given for the problem \
                         is invalid".format(model.ndof))
    
    return ke
