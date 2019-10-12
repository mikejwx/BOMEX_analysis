import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

"""
Plot the profiles of q_v for the BOMEX simlations at different grid spacing.
"""
qv_key = u'STASH_m01s00i010'
z_key  = 'thlev_zsea_theta'

qv_all = {}
exps   = ['DX050', 'Control', 'DX200', 'DX400', 'DX800']
# Read the data
for exp in exps:
    qv_path = '/nerc/n02/n02/xb899100/BOMEX/' + exp + '/water.nc' # path to netCDF containing qv
    qv_nc   = Dataset(qv_path, 'r')
    qv_data = qv_nc.variables[qv_key][:]*1.
    z       = qv_nc.variables[z_key][:]*1.
    qv_nc.close()
    
    qv_all[exp] = np.nanmean(qv_data, axis = (0, 2, 3))

import matplotlib
my_cmap = matplotlib.cm.get_cmap('Reds')
my_colors = {}

for exp in exps:
    my_colors[exp] = [my_cmap((1000-float(exp[2:]))/1000.0) if exp != 'Control' else my_cmap(0.9)][0]

# Initialise the figure
fig = plt.figure()
ax  = fig.add_subplot(1, 1, 1)

for exp in exps:
    ax.plot(qv_all[exp]*1000.0, z/1000.0, color = my_colors[exp], label = [exp if exp != 'Control' else 'REF'][0], lw = 2)

ax.plot([17.0, 16.3, 10.7, 4.2, 3.0], np.array([0., 520., 1480., 2000., 3000.])/1000.0, 'k--', lw = 2, label = 'Initial')
plt.legend(loc = 0)
ax.set_xlim([6, 18])
ax.set_ylim([0, 2.5])
ax.set_xlabel(u'q$_{v}$ (g kg$^{-1}$)')
ax.set_ylabel(u' Height (km)')
plt.savefig('../BOMEX_figs/qv_resolution_comparison.png', dpi = 150)
plt.show()

