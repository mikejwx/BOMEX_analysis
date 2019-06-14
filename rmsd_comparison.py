"""
Code to compare the initial and 'final' conditions between the BOMEX simulations
with different resolutions and domain structures.

Plots a figure with dx along the x-axis and the root mean squared difference 
along the y-axis.

Can compare the 75L and the 140L experiments all in one place.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from BOMEX_tools import *
from scipy import integrate
from multiprocessing import Pool

# Thermodynamic constants (global scope)
Lv = 2.501e6
cpd = 1005.

# List of experiments
exps = ['Control', 'DX050', 'DX050_128', 'DX200', 'DX400', 'DX800', 'L140', 'L140_sp']

# Define the plotting colors
exp_cols = {'DX050'     : 'black',
            'DX050_128' : 'black',
            'DX200'     : 'brown',
            'DX400'     : 'red',
            'DX800'     : 'orange',
            'L140'      : 'b',
            'L140_sp'   : 'b'}

# Define the plotting markers
exp_marker = {'DX050'     : 'o',
              'DX050_128' : '^',
              'DX200'     : 'o',
              'DX400'     : 'o',
              'DX800'     : 'o',
              'L140'      : 'D',
              'L140_sp'   : 's'}

# Assign each experiment to the correct horizontal grid spacing
dx = {'DX050'     : 50.0,
      'DX050_128' : 50.0,
      'Control'   : 100.0,
      'DX200'     : 200.0,
      'DX400'     : 400.0,
      'DX800'     : 800.0,
      'L140'      : 100.0,
      'L140_sp'   : 100.0}

mean_diff_theta = {}
mean_diff_qv    = {}

verbose = False

# Define the initial conditions
init_theta   = np.array([298.7, 298.7, 302.4, 308.2, 311.85])
init_theta_z = np.array([0., 520., 1480., 2000., 3000.])

init_qv   = np.array([17.0, 16.3, 10.7, 4.2, 3.0])
init_qv_z = np.array([0., 520., 1480., 2000., 3000.])

# Read in the data, compute the rmsd, and plot the point for each panel
def do_diff(exp):
    print exp
    if exp in ['Control', 'DX050', 'DX050_128', 'DX200', 'DX400', 'DX800']:
        theta_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/thermo.nc' # path to netCDF containing theta
        qv_path    = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/water.nc' # path to netCDF containing qv
        i_max = None
    elif exp in ['L140', 'L140_sp']:
        theta_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/bouy.nc' # path to netCDF containing theta
        qv_path    = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/mr.nc' # path to netCDF containing qv
    
    # Define some keys to read in variables
    theta_key = u'STASH_m01s00i004'
    qv_key    = u'STASH_m01s00i010'
    
    # Read variables
    theta_nc   = Dataset(theta_path, 'r')
    theta_data = theta_nc.variables[theta_key][:]*1.
    time_key   = [tkey for tkey in theta_nc.variables.keys() if 'min' in tkey][0]
    times      = theta_nc.variables[time_key][:]*1.
    theta_z    = theta_nc.variables['thlev_zsea_theta'][:]*1.
    if exp in ['L140', 'L140_sp']:
        i_max = np.where(np.abs(theta_z - 3000.0) == np.min(np.abs(theta_z - 3000.0)))[0][0]
        if theta_z[i_max] < 3000.0:
            i_max += 1
    
    theta_nc.close()
    
    qv_nc   = Dataset(qv_path, 'r')
    qv_data = qv_nc.variables[qv_key][:]*1.
    qv_nc.close()
    
    # Filter to only hour five of the simulation
    hour5_indexes = [i for i in xrange(len(times)) if 300. <= times[i] <= 360.]
    
    # Compute the mean absolute fractional difference
    init_theta0 = interpolate.interp1d(x = init_theta_z, y = init_theta, fill_value = 'extrapolate')(theta_z[:i_max])
    init_qv0    = interpolate.interp1d(x = init_qv_z,y = init_qv, fill_value = 'extrapolate')(theta_z[:i_max])
    diff_theta = (np.nanmean(theta_data[hour5_indexes,:i_max,:,:], axis = (0, 2, 3)) - init_theta0)/init_theta0
    diff_qv    = (1000*np.nanmean(qv_data[hour5_indexes,:i_max,:,:], axis = (0, 2, 3)) - init_qv0)/init_qv0
    
    mean_diff_theta = integrate.trapz(x = theta_z[:i_max], y = np.abs(diff_theta[:i_max]))/np.max(theta_z[:i_max])
    mean_diff_qv    = integrate.trapz(x = theta_z[:i_max], y = np.abs(diff_qv[:i_max]))/np.max(theta_z[:i_max])
    
    return mean_diff_theta, mean_diff_qv, exp

p = Pool(processes = len(exps))
output = p.map(do_diff, exps)
p.close()
p.join()

# Reconfigure output into the right dictionaries
for out in output:
    mean_diff_theta[out[-1]] = out[0]
    mean_diff_qv[out[-1]] = out[1]

# Start the figure
fs = 14
fig = plt.figure(tight_layout = True, figsize = (14, 8))

# Panel (a) the potential temperature rmsd
ax0 = fig.add_subplot(1, 2, 1)

# Panel (b) the specific humidity rmsd
ax1 = fig.add_subplot(1, 2, 2)

for exp in exps:
    # Add point to plot
    if exp == 'Control':
        ax0.plot(dx[exp], 100*mean_diff_theta[exp], marker = 'o', color = 'b', label = 'REF', ls = 'None', markersize = fs-2)
        ax1.plot(dx[exp], 100*mean_diff_qv[exp], marker = 'o', color = 'b', label = 'REF', ls = 'None', markersize = fs-2)
    else:
        ax0.plot(dx[exp], 100*mean_diff_theta[exp], marker = exp_marker[exp], color = exp_cols[exp], label = exp, ls = 'None', markersize = fs-2)
        ax1.plot(dx[exp], 100*mean_diff_qv[exp], marker = exp_marker[exp], color = exp_cols[exp], label = exp, ls = 'None', markersize = fs-2)
    print 'Finished ' + exp

ax0.legend(loc = 'lower right', ncol = 2, numpoints = 1)
ax0.set_xlabel('Horizontal Grid Spacing (m)', fontsize = fs)
ax0.set_ylabel('$\delta_{0}$, Mean percentage change in $\\theta$', fontsize = fs)
ax0.set_xlim([0, 1000])
ax0.text(50, 0.875*0.3, 'a)', fontsize = fs)

ax1.set_xlabel('Horizontal Grid Spacing (m)', fontsize = fs)
ax1.set_ylabel('$\delta_{0}$, Mean percentage change in q$_{v}$', fontsize = fs)
ax1.set_xlim([0, 1000])
ax1.text(50, 0.875*14, 'b)', fontsize = fs)

plt.savefig('../BOMEX_figs/diff_comparison.png', dpi = 150)
plt.show()

