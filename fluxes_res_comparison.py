"""
Plot the turbulent flux profiles of total water and liquid potential temperature
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
fig = plt.figure(tight_layout = True, figsize = (9, 6))
for exp in exps:
    l_subgrid = True
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
    w_key     = u'STASH_m01s00i150'
    rho_key   = u'STASH_m01s00i389'
    theta_key = u'STASH_m01s00i004'
    qv_key    = u'STASH_m01s00i010'
    mcl_key   = u'STASH_m01s00i392'
    if l_subgrid:
        shf_key   = u'STASH_m01s03i216'
        lhf_key   = u'STASH_m01s03i222'
    
    if verbose:
        print 'Opening netCDFs'
    w_nc   = Dataset(w_path, 'r')
    rho_nc = Dataset(rho_path, 'r')
    if l_subgrid:
        tinc_nc = Dataset(heat_path, 'r')
        qinc_nc = Dataset(moist_path, 'r')
    
    # Get the dimensions for my winds and rho
    if verbose:
        print 'Defining keys for the data dimensions'
    z_w = 'thlev_zsea_theta'
    z_rho = 'rholev_zsea_rho'
    lat_w = 'latitude_t'
    lon_w = 'longitude_t'
    lat_rho = 'latitude_t'
    lon_rho = 'longitude_t'
    
    w_time_shape = w_nc.variables[w_key][:].shape[0]
    w_time_key = [key for key in w_nc.variables.keys() if ('min' in key) and (w_nc.variables[key].shape[0] == w_time_shape)][0]
    if l_subgrid:
        f_time_shape = qinc_nc.variables[lhf_key][:].shape[0]
        f_time_key = [key for key in qinc_nc.variables.keys() if ('min' in key) and (qinc_nc.variables[key].shape[0] == f_time_shape)][0]
    
    if verbose:
        print 'regridding to "rho_data"'
    rho_data = regrid_array(rho_nc.variables[rho_key][:]*1., rho_nc.variables[z_rho][:], rho_nc.variables[lat_rho][:]*1., rho_nc.variables[lon_rho][:]*1., w_nc.variables[z_w][:]*1., w_nc.variables[lat_w][:]*1., w_nc.variables[lon_w][:]*1., verbose)
    
    # Read additional variables
    theta_nc = Dataset(theta_path, 'r')
    theta_data = theta_nc.variables[theta_key][:]*1.
    theta_nc.close()
    
    qv_nc = Dataset(qv_path, 'r')
    qv_data = qv_nc.variables[qv_key][:]*1.
    qv_nc.close()
    
    mcl_nc = Dataset(mcl_path, 'r')
    mcl_data = mcl_nc.variables[mcl_key][:]*1.
    mcl_nc.close()
    
    # Convert to specific humidity form
    qcl_data = mcl_data/(1.0 + mcl_data)

    ### Figure 4 ###
    # Filter to 180 to 360 minutes
    hour3to6_indexes = [i for i in xrange(len(w_nc.variables[w_time_key])) if 180. <= w_nc.variables[w_time_key][i] <= 360.]
    if l_subgrid:
        hr3to6_indexes_f = [i for i in xrange(len(qinc_nc.variables[f_time_key])) if 180. <= qinc_nc.variables[f_time_key][i] <= 360.]
    
    # Put everything on the same vertical levels
    # The subgrid rho levels (z_rho_sg[0] = 0.0m not 20.0m)
    z_rho_sg = np.concatenate((np.array([0.0]), rho_nc.variables[z_rho][1:]))
    z_rho_SG = np.concatenate((np.array([0.0]), rho_nc.variables[z_rho][:])) # keep the level at 20m in there
    lat_rho  = rho_nc.variables[lat_rho][:]*1.
    lon_rho  = rho_nc.variables[lon_rho][:]*1.
    lat_w    = w_nc.variables[lat_w][:]*1.
    lon_w    = w_nc.variables[lon_w][:]*1.
    z_w_data = w_nc.variables[z_w][:]*1.
    
    # Regrid everything
    theta_data_f = regrid_array(theta_data[hour3to6_indexes,:,:,:], z_w_data, lat_w, lon_w, z_rho_SG, lat_rho, lon_rho, verbose)
    qv_data_f    = regrid_array(qv_data[hour3to6_indexes,:,:,:], z_w_data, lat_w, lon_w, z_rho_SG, lat_rho, lon_rho, verbose)
    qcl_data_f   = regrid_array(qcl_data[hour3to6_indexes,:,:,:], z_w_data, lat_w, lon_w, z_rho_SG, lat_rho, lon_rho, verbose)
    w_data_f     = regrid_array(np.concatenate((np.zeros((len(hour3to6_indexes),1,lat_w.shape[0], lon_w.shape[1])), w_nc.variables[w_key][hour3to6_indexes,:,:,:]), axis = 1), np.concatenate((np.array([0.0]), z_w_data)), lat_w, lon_w, z_rho_SG, lat_rho, lon_rho, verbose)
    rho_data_f   = regrid_array(rho_data[hour3to6_indexes,:,:,:], z_w_data, lat_w, lon_w, z_rho_SG, lat_rho, lon_rho, verbose)
    
    if l_subgrid:
        # Read the subgrid theta_l and q_T turbulent fluxes
        shf_data_f = regrid_array(tinc_nc.variables[shf_key][hr3to6_indexes_f,:,:,:], z_rho_sg, lat_rho, lon_rho, z_rho_SG, lat_rho, lon_rho, verbose)
        lhf_data_f = regrid_array(qinc_nc.variables[lhf_key][hr3to6_indexes_f,:,:,:], z_rho_sg, lat_rho, lon_rho, z_rho_SG, lat_rho, lon_rho, verbose)
        tinc_nc.close()
        qinc_nc.close()
    else:
        shf_data_f = np.zeros_like(rho_data_f)
        lhf_data_f = np.zeros_like(rho_data_f)
    
    ### Plot the figure ###
    """
    Figure of the horizontally averaged turbulent fluxes.
    (a) turbulent vertical flux of total water
    (b) turbulent vertical flux of liquid potential temperature
    """
    # Calculate variables
    qt = qv_data_f + qcl_data_f # total water is vapour plus cloud liquid
    qt = rho_data_f*qt/(1. - qt)
    tl = theta_data_f - (Lv/cpd)*qcl_data_f # liquid potential temperature approximation from AMS Glossary
    
    # Calculate the turbulent fluxes
    wp = np.zeros_like(w_data_f)
    wbar = np.mean(w_data_f, axis = (2, 3))
    qtp = np.zeros_like(qt)
    qtbar = np.mean(qt, axis = (2, 3))
    tlp = np.zeros_like(tl)
    tlbar = np.mean(tl, axis = (2, 3))
    
    for it in xrange(w_data_f.shape[0]):
        wp[it,:,:,:]  = np.transpose(np.transpose(w_data_f[it,:,:,:]) - wbar[it,:])
        qtp[it,:,:,:] = np.transpose(np.transpose(qt[it,:,:,:]) - qtbar[it,:])
        tlp[it,:,:,:] = np.transpose(np.transpose(tl[it,:,:,:]) - tlbar[it,:])
    
    wpqtp = np.mean(wp*qtp, axis = (0, 2, 3))
    wptlp = np.mean(rho_data_f*wp*tlp, axis = (0, 2, 3))
    shf_m = np.mean(shf_data_f, axis = (0, 2, 3))
    lhf_m = np.mean(lhf_data_f, axis = (0, 2, 3))
    
    axa = fig.add_subplot(1, 2, 1)
    if exp == 'Control':
        axa.plot(Lv*(wpqtp + lhf_m), z_rho_SG, 'b', lw = 3, label = 'REF')
    else:
        axa.plot(Lv*(wpqtp + lhf_m), z_rho_SG, color = exp_cols[exp], ls = exp_ls[exp], lw = 2, label = exp)
    axa.set_xlim([0, 200])
    axa.set_ylim([0, 2500])
    axa.set_xticks(range(0, 200, 25))
    axa.set_yticks(range(0, 2500, 500))
    axa.set_xlabel("$\overline{w^{\prime}q_{t}^{\prime}}$ (W/m$^{2}$)")
    axa.set_ylabel('height (m)')
    plt.grid()
    
    axb = fig.add_subplot(1, 2, 2)
    if exp == 'Control':
        axb.plot(cpd*wptlp + shf_m, z_rho_SG, 'b', lw = 3, label = 'REF')
    else:
        axb.plot(cpd*wptlp + shf_m, z_rho_SG, color = exp_cols[exp], ls = exp_ls[exp], lw = 2, label = exp)
    axb.set_xlim([-40, 10])
    axb.set_ylim([0, 2500])
    axb.set_xticks(range(-40, 11, 10))
    axb.set_yticks(range(0, 2500, 500))
    axb.set_xlabel("$\overline{w^{\prime}\\theta_{l}^{\prime}}$ (W/m$^{2}$)")
    axb.set_ylabel('height (m)')
    plt.grid()
    
    w_nc.close()
    rho_nc.close()

axa.legend(loc = 0)
plt.savefig('../BOMEX_figs/fluxes_res_comparison.png', dpi = 150)
plt.show()


    
