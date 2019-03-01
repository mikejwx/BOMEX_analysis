import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate

"""
Functions to regrid winds and compute the resolved TKE for our simulations
"""

def regrid_array(in_data, in_z, in_y, in_x, out_z, out_y, out_x, verbose = True):
    """
    Takes a 3- or 4D array of input data and its coordinate axes, then regrids 
    that data to provided coordinate axes. Regridding is performed assuming no 
    change in resolution occurs.
    """
    if verbose:
        print 'Checking which dimensions need to be interpolated'
    # Test 1: Are they the same length?
    z_test = (len(in_z) == len(out_z)) # True = they're the same, False = they're different
    y_test = (len(in_y) == len(out_y))
    x_test = (len(in_x) == len(out_x))
    # Test 2: if they are the same length, do they have the same points?
    if z_test:
        z_test = np.min([in_z[i] == out_z[i] for i in range(len(in_z))])
    if y_test:
        y_test = np.min([in_y[i] == out_y[i] for i in range(len(in_y))])
    if x_test:
        x_test = np.min([in_x[i] == out_x[i] for i in range(len(in_x))])
    
    if len(in_data.shape) == 3:
        if verbose:
            print '"in_data" is a 3D array'
        # Working with in_data[z, y, x]
        # Interpolate in the z-direction first
        if z_test:
            if verbose:
                print 'z-dimension does not need to be interpolated'
            out_data = in_data*1.
        else:
            if verbose:
                print 'Interpolating z-dimension'
            out_data = interpolate.interp1d(x = in_z, y = in_data, axis = 0, fill_value = 'extrapolate')(out_z)
        # Interpolate in the y-direction
        if not y_test:
            if verbose:
                print 'Interpolating y-dimension'
            for i in xrange(in_y.shape[1]):
                 out_data[:,:,i] = interpolate.interp1d(x = in_y[:,i], y = out_data[:,:,i], axis = 1, fill_value = 'extrapolate')(out_y[:,i])
        else:
            if verbose:
                print 'y-dimension does not need to be interpolated'
        # Interpolate in the x-direction
        if not x_test:
            if verbose:
                print 'Interpolating x-dimension'
            for j in xrange(in_x.shape[0]):
                out_data[:,j,:] = interpolate.interp1d(x = in_x[j,:], y = out_data[:,j,:], axis = 1, fill_value = 'extrapolate')(out_x[j,:])
        else:
            if verbose:
                print 'x-dimension does not need to be interpolated'
    elif len(in_data.shape) == 4:
        if verbose:
            print '"in_data" is a 4D array'
        # Working with in_data[t, z, y, x]
        # Interpolate in the z-direction first
        if z_test:
            if verbose:
                print 'z-dimension does not need to be interpolated'
            out_data = in_data*1.
        else:
            if verbose:
                print 'Interpolating the z-dimension'
            out_data = interpolate.interp1d(x = in_z, y = in_data, axis = 1, fill_value = 'extrapolate')(out_z)
        # Interpolate in the y-direction
        if not y_test:
            if verbose:
                print 'Interpolating the y-dimension'
            for i in xrange(in_y.shape[1]):
                out_data[:,:,:,i] = interpolate.interp1d(x = in_y[:,i], y = out_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(out_y[:,i])
        else:
            if verbose:
                print 'y-dimension does not need to be interpolated'
        # Interpolate in the x-direction
        if not x_test:
            if verbose:
                print 'Interpolating the x-dimension'
            for j in xrange(in_x.shape[0]):
                out_data[:,:,j,:] = interpolate.interp1d(x = in_x[j,:], y = out_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(out_x[j,:])
        else:
            if verbose:
                print 'x-dimension does not need to be interpolated'
    
    return out_data




    # 3. interpolate in y
    if target_y_key != current_y_key:
        print 'Regridding in latitude...'
        # Interpolate in y
        for i in xrange(current_data.shape[3]):
            current_data[:,:,:,i] = interpolate.interp1d(current_y[:,i], current_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(target_y[:,i])
        print 'Complete.'
    # 4. interpolate in x
    if target_x_key != current_x_key:
        print 'Regridding in longitude...'
        # Interpolate in x
        for j in xrange(current_data.shape[2]):
            current_data[:,:,j,:] = interpolate.interp1d(current_x[j,:], current_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(target_x[j,:])
        print 'Complete.'

def get_tke(u, v, w, rho = 1., verbose = True):
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
        if verbose:
            print 'This is a 2D calculation'
        # 2D array
        TKE = 0.5*rho*((u - np.mean(u))**2. + (v - np.mean(v))**2. + (w - np.mean(w))**2.)
    elif nd == 3:
        if verbose:
            print 'This is a 3D calculation'
        # 3D array
        up = np.transpose(np.transpose(u) - np.mean(u, axis = (1, 2)))
        vp = np.transpose(np.transpose(v) - np.mean(v, axis = (1, 2)))
        wp = np.transpose(np.transpose(w) - np.mean(w, axis = (1, 2)))
        TKE = 0.5*rho*(up**2. + vp**2. + wp**2.)
    elif nd == 4:
        if verbose:
            print 'This is a 4D calculation'
            print 'Calculating the wind perturbations'
        # 4D array
        up = np.zeros_like(u)
        vp = np.zeros_like(v)
        wp = np.zeros_like(w)
        u_mean = np.mean(u, axis = (2, 3))
        v_mean = np.mean(v, axis = (2, 3))
        w_mean = np.mean(w, axis = (2, 3))
        for j in xrange(u.shape[2]):
            for i in xrange(u.shape[3]):
                if verbose:
                    print '[j,i] = [' + str(j) + ',' + str(i) + ']'
                up[:,:,j,i] = u[:,:,j,i] - u_mean
                vp[:,:,j,i] = v[:,:,j,i] - v_mean
                wp[:,:,j,i] = w[:,:,j,i] - w_mean
        TKE = 0.5*rho*(up**2. + vp**2. + wp**2.)
    return TKE
