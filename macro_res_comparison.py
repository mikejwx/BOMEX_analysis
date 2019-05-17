"""
Plot the total cloud cover, liquid water path, and turbulent kinetic energy
for each of the resolutions of the L75 BOMEX experiments.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from BOMEX_tools import *
from scipy import integrate

# Thermodynamic constants (global scope)
Lv = 2.501e6
cpd = 1005.

# List of experiments
exps = ['Control', 'DX050', 'DX050_128', 'DX200', 'DX400', 'DX800']
exp_cols = {'DX050' : 'black',
            'DX050_128' : 'black',
            'DX100' : 'brown',
            'DX200' : 'red',
            'DX400' : 'orange',
            'DX800' : 'pink'}
exp_ls   = {'DX050' : '-',
            'DX050_128' : '--',
            'DX100' : '-',
            'DX200' : '-',
            'DX400' : '-',
            'DX800' : '-'}
verbose = False
fig = plt.figure(tight_layout = True, figsize = (6, 8))
ax0 = fig.add_subplot(3, 1, 1)
ax0.set_xlim([0, 360])
a_max = 0.3
a_inc = 0.1
ax0.set_ylim([0, a_max])
ax0.set_xticks(np.arange(0, 361, 60))
ax0.set_yticks(np.arange(0, a_max + 0.01, a_inc))
ax0.set_ylabel('total cloud cover')
ax0.text(10, 0.875*a_max, 'a)')

ax1 = fig.add_subplot(3, 1, 2)
ax1.set_xlim([0, 360])
b_max = 100
b_inc = 25.
ax1.set_ylim([0, b_max])
ax1.set_xticks(np.arange(0, 361, 60))
ax1.set_yticks(np.arange(0, b_max + 1, b_inc))
ax1.set_ylabel('LWP (gm$^{-2}$)')
ax1.text(10, 0.875*b_max, 'b)')

ax2 = fig.add_subplot(3, 1, 3)
ax2.set_xlim([0, 360])
c_max = 1250
c_inc = 250.
ax2.set_ylim([0, c_max])
ax2.set_xticks(np.arange(0, 361, 60))
ax2.set_yticks(np.arange(0, c_max+1, c_inc))
ax2.set_xlabel('Time (min)')
ax2.set_ylabel('TKE (kgm$^{-1}$s$^{-2}$)')
ax2.text(10, 625, 'c)')

for exp in exps:
    lwp_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/water.nc' # path to netCDF containing LWP for L75 sims
    u_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/wind.nc' # path to netCDF containing u
    v_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/wind.nc' # path to netCDF containing v
    w_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/wind.nc' # path to netCDF containing w
    rho_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/thermo.nc' # path to netCDF containing rho
    theta_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/thermo.nc' # path to netCDF containing theta
    qv_path    = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/water.nc' # path to netCDF containing qv
    mcl_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/water.nc' # path to netCDF containing mcl
    heat_path  = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/tempinc.nc' # path to netCDF containing boundary layer heat fluxes (subgrid?)
    moist_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/qinc.nc' # path to netCDF containing boundary layer moisture fluxes (subgrid?)
    
    
    # Read in the netCDF data and process it to create variables for our plots
    if verbose:
        print 'Defining keys for data variables'
    lwp_key   = u'STASH_m01s30i405'
    u_key     = u'STASH_m01s00i002'
    v_key     = u'STASH_m01s00i003'
    w_key     = u'STASH_m01s00i150'
    rho_key   = u'STASH_m01s00i389'
    
    if verbose:
        print 'Opening netCDFs'
    lwp_nc = Dataset(lwp_path, 'r')
    u_nc   = Dataset(u_path, 'r')
    v_nc   = Dataset(v_path, 'r')
    w_nc   = Dataset(w_path, 'r')
    rho_nc = Dataset(rho_path, 'r')
    
    # Get the dimensions for my winds and rho
    if verbose:
        print 'Defining keys for the data dimensions'
    z_u = 'rholev_zsea_rho'
    z_v = 'rholev_zsea_rho'
    z_w = 'thlev_zsea_theta'
    z_rho = 'rholev_zsea_rho'
    lat_u = 'latitude_cu'
    lon_u = 'longitude_cu'
    lat_v = 'latitude_cv'
    lon_v = 'longitude_cv'
    lat_w = 'latitude_t'
    lon_w = 'longitude_t'
    lat_rho = 'latitude_t'
    lon_rho = 'longitude_t'
    
    w_time_shape = w_nc.variables[w_key][:].shape[0]
    w_time_key = [key for key in w_nc.variables.keys() if ('min' in key) and (w_nc.variables[key].shape[0] == w_time_shape)][0]
    
    if verbose:
        print 'Regridding the data for TKE calculation'
        print 'regridding to "u_data"'
    u_data = regrid_array(u_nc.variables[u_key][:]*1., u_nc.variables[z_u][:]*1., u_nc.variables[lat_u][:]*1., u_nc.variables[lon_u][:]*1., w_nc.variables[z_w][:], w_nc.variables[lat_w][:], w_nc.variables[lon_w][:]*1., verbose)
    if verbose:
        print 'regridding to "v_data"'
    v_data = regrid_array(v_nc.variables[v_key][:]*1., v_nc.variables[z_v][:]*1., v_nc.variables[lat_v][:]*1., v_nc.variables[lon_v][:]*1., w_nc.variables[z_w][:], w_nc.variables[lat_w][:], w_nc.variables[lon_w][:]*1., verbose)
    if verbose:
        print 'regridding to "rho_data"'
    rho_data = regrid_array(rho_nc.variables[rho_key][:]*1., rho_nc.variables[z_rho][:], rho_nc.variables[lat_rho][:]*1., rho_nc.variables[lon_rho][:]*1., w_nc.variables[z_w][:]*1., w_nc.variables[lat_w][:]*1., w_nc.variables[lon_w][:]*1., verbose)
    
    if verbose:
        print 'Calculating the TKE'
    TKE_time = w_nc.variables[w_time_key][:]
    TKE = get_tke(u_data, v_data, w_nc.variables[w_key][:]*1., rho_data, verbose)
    # Compute the vertical integral
    if verbose:
        print 'Integrating the TKE'
    viTKE = np.mean(integrate.trapz(y = TKE, x = w_nc.variables[z_w][:]*1., axis = 1), axis = (1,2))
    lwp_time_shape = lwp_nc.variables[lwp_key][:].shape[0]
    lwp_time_key = [key for key in lwp_nc.variables.keys() if ('min' in key) and (lwp_nc.variables[key].shape[0] == lwp_time_shape)][0]
    
    # Panel (a)
    if exp == 'Control':
        ax0.plot(lwp_nc.variables[lwp_time_key][:]*1., np.nanmean(np.where((lwp_nc.variables[lwp_key][:]*1. >= 5./1000.), 1., 0.), axis = (1, 2)), color = 'b', lw = 3, label = 'REF')
    else:
        ax0.plot(lwp_nc.variables[lwp_time_key][:]*1., np.nanmean(np.where((lwp_nc.variables[lwp_key][:]*1. >= 5./1000.), 1., 0.), axis = (1, 2)), color = exp_cols[exp], ls = exp_ls[exp], lw = 2, label = exp)
    # Panel (b)
    if exp == 'Control':
        ax1.plot(lwp_nc.variables[lwp_time_key][:]*1., 1000.*np.nanmean(lwp_nc.variables[lwp_key][:]*1., axis = (1, 2)), color = 'b', lw = 3, label = 'REF')
    else:
        ax1.plot(lwp_nc.variables[lwp_time_key][:]*1., 1000.*np.nanmean(lwp_nc.variables[lwp_key][:]*1., axis = (1, 2)), color = exp_cols[exp], ls = exp_ls[exp], lw = 2, label = exp)
    # Panel (c)
    if exp == 'Control':
        ax2.plot(TKE_time, viTKE, color = 'b', lw = 3, label = 'REF')
    else:
        ax2.plot(TKE_time, viTKE, color = exp_cols[exp], ls = exp_ls[exp], lw = 2, label = exp)
    
    w_nc.close()
    rho_nc.close()

ax0.legend(loc = 0)
plt.savefig('../BOMEX_figs/macro_res_comparison.png', dpi = 150)
plt.show()


