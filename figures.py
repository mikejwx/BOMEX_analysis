import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from BOMEX_tools import *
from scipy import integrate

"""
Defining functions to plot the figures from Siebesma et al. (2003)
Figure 2: Time series of (a) the total cloud cover, (b) the vertically 
          integrated cloud liquid water path, and (c) the vertically integrated 
          TKE.

Figure 3: Mean profiles averaged over the last hour of (a) potential 
          temperature, theta, (b) water vapour specific humidity q_v, (c) the 
          horizontal velocity components, and (d) the liquid water q_l.

Figure 4: Turbulent flux profiles averaged over the last 3 h of (a) total water 
          q_t, (b) liquid water potential temperature theta_l, (c) liquid water 
          q_l, (d) virtual potential temperature theta_v, and (e) zonal wind u.

Figure 5: Profiles of (a) turbulent kinetic energy (TKE) and its vertical 
          component (b) sigma_w^2.
"""

def figure2(lwp_nc, lwp_key, lwp_time_key, viTKE, TKE_time, fig_name):
    """
    Compute the time series and plot them.
    """
    fig = plt.figure(figsize=(6, 8), tight_layout = True)
    # Panel (a)
    ax0 = fig.add_subplot(3, 1, 1)
    ax0.plot(lwp_nc.variables[lwp_time_key][:]*1., np.nanmean(np.where((lwp_nc.variables[lwp_key][:]*1. >= 5./1000.), 1., 0.), axis = (1, 2)), color = 'k', lw = 2)
    ax0.set_xlim([0, 360])
    ax0.set_ylim([0, 0.2])
    ax0.set_xticks(np.arange(0, 361, 60))
    ax0.set_yticks(np.arange(0, 0.21, 0.1))
    ax0.set_ylabel('total cloud cover')
    ax0.text(10, 0.175, 'a)')
    # Panel (b)
    ax1 = fig.add_subplot(3, 1, 2)
    ax1.plot(lwp_nc.variables[lwp_time_key][:]*1., 1000.*np.nanmean(lwp_nc.variables[lwp_key][:]*1., axis = (1, 2)), color = 'k', lw = 2)
    ax1.set_xlim([0, 360])
    ax1.set_ylim([0, 20])
    ax1.set_xticks(np.arange(0, 361, 60))
    ax1.set_yticks(np.arange(0, 21, 5))
    ax1.set_ylabel('LWP (gm$^{-2}$)')
    ax1.text(10, 17.5, 'b)')
    # Panel (c)
    ax2 = fig.add_subplot(3, 1, 3)
    ax2.plot(TKE_time, viTKE, color = 'k', lw = 2)
    ax2.set_xlim([0, 360])
    ax2.set_ylim([0, 750])
    ax2.set_xticks(np.arange(0, 361, 60))
    ax2.set_yticks(np.arange(0, 751, 250))
    ax2.set_xlabel('Time (min)')
    ax2.set_ylabel('TKE (kgm$^{-1}$s$^{-2}$)')
    ax2.text(10, 625, 'c)')
    
    plt.savefig('../Figure2_Siebesma2003_'+ fig_name + '.png', dpi = 150)
    plt.show()

def main(lwp_path, u_path, v_path, w_path, rho_path, fig_name, verbose = True):
    """
    Read in the netCDF data and process it to create variables for our plots
    """
    if verbose:
        print 'Defining keys for data variables'
    lwp_key = u'STASH_m01s30i405'
    u_key   = u'STASH_m01s00i002'
    v_key   = u'STASH_m01s00i003'
    w_key   = u'STASH_m01s00i150'
    rho_key = u'STASH_m01s00i389'
    
    if verbose:
        print 'Opening netCDFs'
    lwp_nc = Dataset(lwp_path, 'r')
    u_nc = Dataset(u_path, 'r')
    v_nc = Dataset(v_path, 'r')
    w_nc = Dataset(w_path, 'r')
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
    
    if verbose:
        print 'Plotting Figure 2'
    figure2(lwp_nc, lwp_key, lwp_time_key, viTKE, TKE_time, fig_name)
    lwp_nc.close()
    u_nc.close()
    v_nc.close()
    w_nc.close()
    rho_nc.close()
    
lwp_path = '/nerc/n02/n02/xb899100/BOMEX/Control/water.nc'
u_path = '/nerc/n02/n02/xb899100/BOMEX/Control/wind.nc'
v_path = '/nerc/n02/n02/xb899100/BOMEX/Control/wind.nc'
w_path = '/nerc/n02/n02/xb899100/BOMEX/Control/wind.nc'
rho_path = '/nerc/n02/n02/xb899100/BOMEX/Control/thermo.nc'
my_fig_name = 'Control'
main(lwp_path, u_path, v_path, w_path, rho_path, my_fig_name, verbose = False)

