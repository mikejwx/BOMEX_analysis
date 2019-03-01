import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate

"""
Functions to regrid winds and compute the resolved TKE for our simulations
"""

def regrid_array(in_data, in_z, in_y, in_x, out_z, out_y, out_x):
    """
    Takes a 3- or 4D array of input data and its coordinate axes, then regrids 
    that data to provided coordinate axes. Regridding is performed assuming no 
    change in resolution occurs.
    """
    # Test 1: Are they the same length?
    z_test = (len(in_z) == len(out_z)) # True = they're the same, False = they're different
    y_test = (len(in_y) == len(out_y))
    x_test = (len(in_x) == len(out_x))
    # Test 2: if they are the same length, do they have the same points?
    if z_test:
        z_test = np.min([in_z[i] == out_z[i] for i in len(in_z)])
    if y_test:
        y_test = np.min([in_y[i] == out_y[i] for i in len(in_y)])
    if x_test:
        x_test = np.min([in_x[i] == out_x[i] for i in len(in_x)])
    
    if len(in_data.shape) == 3:
        # Working with in_data[z, y, x]
        # Interpolate in the z-direction first
        if z_test:
            out_data = in_data*1.
        else:
            out_data = interpolate.interp1d(x = in_z, y = in_data, axis = 0, fill_value = 'extrapolate')(out_z)
        # Interpolate in the y-direction
        if y_test:
            out_data = interpolate.interp1d(x = in_y, y = out_data, axis = 1, fill_value = 'extrapolate')(out_y)
        # Interpolate in the x-direction
        if x_test:
            out_data = interpolate.interp1d(x = in_x, y = out_data, axis = 2, fill_value = 'extrapolate')(out_x)
    elif len(in_data.shape) == 4:
        # Working with in_data[t, z, y, x]
        # Interpolate in the z-direction first
        if z_test:
            out_data = in_data*1.
        else:
            out_data = interpolate.interp1d(x = in_z, y = in_data, axis = 1, fill_value = 'extrapolate')(out_z)
        # Interpolate in the y-direction
        if y_test:
            out_data = interpolate.interp1d(x = in_y, y = out_data, axis = 2, fill_value = 'extrapolate')(out_y)
        # Interpolate in the x-direction
        if x_test:
            out_data = interpolate.interp1d(x = in_x, y = out_data, axis = 3, fill_value = 'extrapolate')(out_x)
    
    return out_data

def get_tke(u, v, w, rho = 1.):
    """
    Function to compute the resolved TKE for u, v, and w wind components given 
    that they are all on the same grid.
    ----------------------------------------------------------------------------
    INPUT:
    u, v, w = the wind component in the x-, y-, and z- directions respectively
    These could be 2-, 3-, or 4D arrays
    Either [y, x], [z, y, x], or [t, z, y, x]
    But they must all be the same shape!
    rho = the air density
    If rho is supplied, it too must be on the same shape/grid as the winds
    OUTPUT:
    TKE = an array the same shape as the input containing the TKE
    
    TKE = (rho/2)*(u'^2 + v'^2 + w'^2) where e.g. u' is a perturbation from
          the horizontal mean
    if rho supplied:
        TKE has units of kg m^2 s^-2
    if rho not supplied:
        TKE has units of m^2 s^-2
    """
    nd = len(u.shape)
    if nd == 2:
        # 2D array
        TKE = 0.5*rho*((u - np.mean(u))**2. + (v - np.mean(v))**2. + (w - np.mean(w))**2.)
    elif nd == 3:
        # 3D array
        TKE = 0.5*rho*((u - np.mean(u, axis = (1, 2)))**2. + (v - np.mean(v, axis = (1, 2)))**2. + (w - np.mean(w, axis = (1, 2)))**2.)
    elif nd == 4:
        # 4D array
        TKE = 0.5*rho*((u - np.mean(u, axis = (2, 3)))**2. + (v - np.mean(v, axis = (2, 3)))**2. + (w - np.mean(w, axis = (2, 3)))**2.)
    return TKE
