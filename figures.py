import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

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

def figure2(lwp_nc, lwp_key, lwp_time_key):
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
    #ax2.plot(tke_variables[time_key][:]*1., np.nanmean(viTKE, axis = (1, 2)), color = 'k', lw = 2)
    ax2.set_xlim([0, 360])
    ax2.set_ylim([0, 750])
    ax2.set_xticks(np.arange(0, 361, 60))
    ax2.set_yticks(np.arange(0, 751, 250))
    ax2.set_xlabel('Time (min)')
    ax2.set_ylabel('TKE (kgm$^{-1}$s$^{-2}$)')
    ax2.text(10, 625, 'c)')
    
    plt.show()
    #plt.savefig('../Figure2_Siebesma2003.png', dpi = 150)

def main(path):
    """
    Read in the netCDF data and process it to create variables for our plots
    """
    lwp_key = u'STASH_m01s30i405'
    lwp_nc = Dataset(path + 'water.nc', 'r')
    lwp_time_shape = lwp_nc.variables[lwp_key][:].shape[0]
    lwp_time_key = [key for key in lwp_nc.variables.keys() if ('min' in key) and (lwp_nc.variables[key].shape[0] == lwp_time_shape)][0]
    figure2(lwp_nc, lwp_key, lwp_time_key)

path = '/nerc/n02/n02/xb899100/BOMEX/Control/'
main(path)
