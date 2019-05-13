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
    
    plt.savefig('../BOMEX_figs/Figure2_Siebesma2003_'+ fig_name + '.png', dpi = 150)
    #plt.show()

def figure3(theta, qv, ql, u, v, z, fig_name):
    """
    Compute the horizonal and time mean profiles for available data from T = 300 
    to T = 360.
    (a) Potential temperature
    (b) Specific Humidity
    (c) u and v wind components
    (d) specific humidity
    Assumes that you pass the variables at the correct simulated times, and that
    all of the data are on the same vertical grid.
    """
    theta_mean = np.nanmean(theta, axis = (0, 2, 3))
    theta_std = np.std(theta, axis = (0, 2, 3))
    qv_mean = np.nanmean(qv, axis = (0, 2, 3))*1000.
    qv_std = np.std(qv, axis = (0, 2, 3))*1000.
    ql_mean = np.nanmean(ql, axis = (0, 2, 3))*1000.
    ql_std = np.std(ql, axis = (0, 2, 3))*1000.
    u_mean = np.nanmean(u, axis = (0, 2, 3))
    u_std = np.std(u, axis = (0, 2, 3))
    v_mean = np.nanmean(v, axis = (0, 2, 3))
    v_std = np.std(v, axis = (0, 2, 3))
    
    fig = plt.figure(figsize = (9, 9), tight_layout = True)
    # Potential Temperature
    ax0 = fig.add_subplot(2, 2, 1)
    ax0.plot([298.7, 298.7, 302.4, 308.2, 311.85], [0., 520., 1480., 2000., 3000.], 'r--')
    ax0.plot(theta_mean, z, color = 'k', lw = 2)
    ax0.fill_betweenx(z, theta_mean - theta_std, theta_mean + theta_std, color = 'k', alpha = 0.5, edgecolor = '')
    ax0.set_xlim([298, 310])
    ax0.set_ylim([0, 2500])
    ax0.text(299, 2250, 'a)')
    ax0.set_xlabel('$\\theta$ (K)')
    ax0.set_ylabel('height (m)')
    # Specific humidity
    ax1 = fig.add_subplot(2, 2, 2)
    ax1.plot([17.0, 16.3, 10.7, 4.2, 3.0], [0., 520., 1480., 2000., 3000.], 'r--')
    ax1.plot(qv_mean, z, color = 'k', lw = 2)
    ax1.fill_betweenx(z, qv_mean - qv_std, qv_mean + qv_std, color = 'k', alpha = 0.5, edgecolor = '')
    ax1.set_xlim([4, 18])
    ax1.set_ylim([0, 2500])
    ax1.text(5, 2250, 'b)')
    ax1.set_xlabel('q$_v$ (g/kg)')
    ax1.set_ylabel('height (m)')
    # Winds
    ax2 = fig.add_subplot(2, 2, 3)
    ax2.plot([-8.75, -8.75, -4.61], [0., 700., 3000.], 'r--')
    ax2.plot([0., 0.], [0., 3000.], 'b--')
    ax2.plot(u_mean, z, color = 'k', lw = 2)
    ax2.fill_betweenx(z, u_mean - u_std, u_mean + u_std, color = 'k', alpha = 0.5, edgecolor = '')
    ax2.plot(v_mean, z, color = 'k', lw = 2)
    ax2.fill_betweenx(z, v_mean - v_std, v_mean + v_std, color = 'k', alpha = 0.5, edgecolor = '')
    ax2.set_xlim([-10, 2])
    ax2.set_ylim([0, 2500])
    ax2.text(-8.75, 2250, 'c)')
    ax2.set_xlabel('u (m/s)    v (m/s)')
    ax2.set_ylabel('height (m)')
    # Cloud liquid
    ax3 = fig.add_subplot(2, 2, 4)
    ax3.plot(ql_mean, z, color = 'k', lw = 2)
    ax3.set_xlim([0, 0.01])
    ax3.set_ylim([0, 2500])
    ax3.text(0.001, 2250, 'b)')
    ax3.set_xlabel('q$_l$ (g/kg)')
    ax3.set_ylabel('height (m)')
    
    plt.savefig('../BOMEX_figs/Figure3_Siebesma2003_' + fig_name + '.png', dpi = 150)
    #plt.show()

def figure4(theta, qv, ql, u, w, rho, shf, lhf, z, fig_name):
    """
    Figure of the horizontally averaged turbulent fluxes.
    (a) turbulent vertical flux of total water
    (b) turbulent vertical flux of liquid potential temperature
    (c) turbulent vertical flux of cloud liquid water
    (d) turbulent vertical flux og virtual potential temperature
    (e) turbulent zonal momentum flux
    ----------------------------------------------------------------------------
    INPUT:
    Requires qv, ql, theta, u, and w on the same grid in 4D f[t,z,y,x]
    only input the correct times, will be mean across whole time axis
    
    will plot just the resolved turbulent fluxes for now.
    """
    Lv = 2.501e6
    cpd = 1005.
    # Calculate variables
    qt = qv + ql # total water is vapour plus cloud liquid
    tl = theta - (Lv/cpd)*ql # liquid potential temperature approximation from AMS Glossary
    tv = theta*(1. + 0.608*qv) # virtual potential temperature
    
    # Calculate the turbulent fluxes
    wp = np.zeros_like(w)
    wbar = np.mean(w, axis = (0, 2, 3))
    qtp = np.zeros_like(qt)
    qtbar = np.mean(qt, axis = (0, 2, 3))
    tlp = np.zeros_like(tl)
    tlbar = np.mean(tl, axis = (0, 2, 3))
    qlp = np.zeros_like(ql)
    qlbar = np.mean(ql, axis = (0, 2, 3))
    tvp = np.zeros_like(tv)
    tvbar = np.mean(tv, axis = (0, 2, 3))
    up = np.zeros_like(u)
    ubar = np.mean(u, axis = (0, 2, 3))
    
    for it in xrange(w.shape[0]):
        wp[it,:,:,:]  = np.transpose(np.transpose(w[it,:,:,:]) - wbar)
        qtp[it,:,:,:] = np.transpose(np.transpose(qt[it,:,:,:]) - qtbar)
        tlp[it,:,:,:] = np.transpose(np.transpose(tl[it,:,:,:]) - tlbar)
        qlp[it,:,:,:] = np.transpose(np.transpose(ql[it,:,:,:]) - qlbar)
        tvp[it,:,:,:] = np.transpose(np.transpose(tv[it,:,:,:]) - tvbar)
        up[it,:,:,:]  = np.transpose(np.transpose(u[it,:,:,:]) - ubar)
    
    wpqtp = np.mean(rho*wp*qtp + lhf, axis = (0, 2, 3))
    wptlp = np.mean(rho*wp*tlp, axis = (0, 2, 3))
    shf_m = np.mean(shf, axis = (0, 2, 3))
    wpqlp = np.mean(rho*wp*qlp, axis = (0, 2, 3))
    wptvp = np.mean(rho*wp*tvp, axis = (0, 2, 3))
    wpup  = np.mean(wp*up, axis = (0, 2, 3))
    
    fig = plt.figure(tight_layout = True, figsize = (9, 12))
    axa = fig.add_subplot(3, 2, 1)
    axa.plot(Lv*wpqtp, z, 'k')
    axa.set_xlim([0, 175])
    axa.set_ylim([0, 2500])
    axa.set_xticks(range(0, 175, 25))
    axa.set_yticks(range(0, 2500, 500))
    axa.set_xlabel("$\overline{w^{\prime}q_{t}^{\prime}}$ (W/m$^{2}$)")
    axa.set_ylabel('height (m)')
    plt.grid()
    
    axb = fig.add_subplot(3, 2, 2)
    axb.plot(cpd*wptlp + shf_m, z, 'k')
    axb.set_xlim([-40, 10])
    axb.set_ylim([0, 2500])
    axb.set_xticks(range(-40, 11, 10))
    axb.set_yticks(range(0, 2500, 500))
    axb.set_xlabel("$\overline{w^{\prime}\\theta_{l}^{\prime}}$ (W/m$^{2}$)")
    axb.set_ylabel('height (m)')
    plt.grid()
    
    axc = fig.add_subplot(3, 2, 3)
    axc.plot(Lv*wpqlp, z, 'k')
    axc.set_xlim([0, 40])
    axc.set_ylim([0, 2500])
    axc.set_xticks(range(0, 41, 10))
    axc.set_yticks(range(0, 2500, 500))
    axc.set_xlabel("$\overline{w^{\prime}q_{l}^{\prime}}$ (W/m$^{2}$)")
    axc.set_ylabel('height (m)')
    plt.grid()
    
    axd = fig.add_subplot(3, 2, 4)
    axd.plot(cpd*wptvp, z, 'k')
    axd.set_xlim([-10, 30])
    axd.set_ylim([0, 2500])
    axd.set_xticks(range(-10, 31, 10))
    axd.set_yticks(range(0, 2500, 500))
    axd.set_xlabel("$\overline{w^{\prime}\\theta_{v}^{\prime}}$ (W/m$^{2}$)")
    axd.set_ylabel('height (m)')
    plt.grid()
    
    axe = fig.add_subplot(3, 2, 5)
    axe.plot(wpup, z, 'k')
    axe.set_xlim([0.05, 0.1])
    axe.set_ylim([0, 2500])
    axe.set_xticks(np.arange(-0.05, 0.11, 0.05))
    axe.set_yticks(range(0, 2500, 500))
    axe.set_xlabel("$\overline{w^{\prime}u^{\prime}}$ (m$^{2}$/s$^{2}$)")
    axe.set_ylabel('height (m)')
    plt.grid()
    
    plt.savefig('../BOMEX_figs/Figure4_Siebesma2003_'+fig_name+'.png', dpi = 150)
    plt.show()

def main(lwp_path, u_path, v_path, w_path, rho_path, theta_path, qv_path, mcl_path, fig_name, verbose = True):
    """
    Read in the netCDF data and process it to create variables for our plots
    """
    if verbose:
        print 'Defining keys for data variables'
    lwp_key   = u'STASH_m01s30i405'
    u_key     = u'STASH_m01s00i002'
    v_key     = u'STASH_m01s00i003'
    w_key     = u'STASH_m01s00i150'
    rho_key   = u'STASH_m01s00i389'
    theta_key = u'STASH_m01s00i004'
    qv_key    = u'STASH_m01s00i010'
    mcl_key   = u'STASH_m01s00i392'
    shf_key   = u'STASH_m01s03i216'
    lhf_key   = u'STASH_m01s03i222'
    
    if verbose:
        print 'Opening netCDFs'
    lwp_nc = Dataset(lwp_path, 'r')
    u_nc   = Dataset(u_path, 'r')
    v_nc   = Dataset(v_path, 'r')
    w_nc   = Dataset(w_path, 'r')
    rho_nc = Dataset(rho_path, 'r')
    tinc_nc = Dataset(heat_path, 'r')
    qinc_nc = Dataset(moist_path, 'r')
    
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
    f_time_shape = qinc_nc.variables[lhf_key][:].shape[0]

    w_time_key = [key for key in w_nc.variables.keys() if ('min' in key) and (w_nc.variables[key].shape[0] == w_time_shape)][0]
    f_time_key = [key for key in qinc_nc.variables.keys() if ('min' in key) and (qinc_nc.variables[key].shape[0] == f_time_shape)][0]
    
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
    
    # read additional variables for figure 3
    theta_nc = Dataset(theta_path, 'r')
    theta_data = theta_nc.variables[theta_key][:]*1.
    theta_nc.close()
    
    qv_nc = Dataset(qv_path, 'r')
    qv_data = qv_nc.variables[qv_key][:]*1.
    qv_nc.close()
    
    mcl_nc = Dataset(mcl_path, 'r')
    mcl_data = mcl_nc.variables[mcl_key][:]*1.
    mcl_nc.close()
    
    # Filter to only hour five of the simulation
    hour5_indexes = [i for i in xrange(len(w_nc.variables[w_time_key])) if 300. <= w_nc.variables[w_time_key][i] <= 360.]
    figure3(theta_data[hour5_indexes,:,:,:], qv_data[hour5_indexes,:,:,:], mcl_data[hour5_indexes,:,:,:], u_data[hour5_indexes,:,:,:], v_data[hour5_indexes,:,:,:], w_nc.variables[z_w][:]*1., fig_name)
    
    # Filter to 180 to 360 minutes
    hour3to6_indexes = [i for i in xrange(len(w_nc.variables[w_time_key])) if 180. <= w_nc.variables[w_time_key][i] <= 360.]
    hr3to6_indexes_f = [i for i in xrange(len(qinc_nc.variables[f_time_key])) if 180. <= qinc_nc.variables[f_time_key][i] <= 360.]
    
    qcl_data = mcl_data[hour3to6_indexes,:,:,:]/(1. + mcl_data[hour3to6_indexes,:,:,:])
    # Read the subgrid theta_l and q_T turbulent fluxes
    shf_data = regrid_array(tinc_nc.variables[shf_key][:]*1., tinc_nc.variables[z_rho][:], tinc_nc.variables[lat_rho][:]*1., tinc_nc.variables[lon_rho][:]*1., w_nc.variables[z_w][:]*1., w_nc.variables[lat_w][:]*1., w_nc.variables[lon_w][:]*1., verbose)
    lhf_data = regrid_array(qinc_nc.variables[lhf_key][:]*1., qinc_nc.variables[z_rho][:], qinc_nc.variables[lat_rho][:]*1., qinc_nc.variables[lon_rho][:]*1., w_nc.variables[z_w][:]*1., w_nc.variables[lat_w][:]*1., w_nc.variables[lon_w][:]*1., verbose)
    tinc_nc.close()
    qinc_nc.close()
    
    figure4(theta_data[hour3to6_indexes,:,:,:], qv_data[hour3to6_indexes,:,:,:], qcl_data, u_data[hour3to6_indexes,:,:,:], w_nc.variables[w_key][hour3to6_indexes,:,:,:]*1., rho_data[hour3to6_indexes,:,:,:], shf_data[hr3to6_indexes_f,:,:,:], lhf_data[hr3to6_indexes_f,:,:,:], w_nc.variables[z_w][:]*1., fig_name)
    lwp_nc.close()
    u_nc.close()
    v_nc.close()
    w_nc.close()
    rho_nc.close()
    
# Define the path to the netCDF containing our variables
exp = 'L140_noBL'
l_control = 0
if l_control:
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
else:
    lwp_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/lwp.nc' # path to netCDF containing LWP L140 sims
    u_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/u.nc' # path to netCDF containing u
    v_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/v.nc' # path to netCDF containing v
    w_path     = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/bouy.nc' # path to netCDF containing w
    rho_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/fluxes.nc' # path to netCDF containing rho
    theta_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/bouy.nc' # path to netCDF containing theta
    qv_path    = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/mr.nc' # path to netCDF containing qv
    mcl_path   = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/mr.nc' # path to netCDF containing mcl
    heat_path  = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/tempinc.nc' # path to netCDF containing boundary layer heat fluxes (subgrid?)
    moist_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/qinc.nc' # path to netCDF containing boundary layer moisture fluxes (subgrid?)

main(lwp_path, u_path, v_path, w_path, rho_path, theta_path, qv_path, mcl_path, exp, verbose = 0)

