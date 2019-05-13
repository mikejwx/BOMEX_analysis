"""
Define and plot the initial conditions for the BOMEX simulation.
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the initial conditions
zthermo = np.array([0, 520, 1480, 2000, 3000])
qt = np.array([17.0, 16.3, 10.7, 4.2, 3.0])
thetal = np.array([298.7, 298.7, 302.4, 308.2, 311.85])

zwind = np.array([0, 700, 3000])
u = np.array([-8.75, -8.75, -4.61])
v = np.array([0, 0, 0])

# Define the large-scale forcing
zgeo = np.array([0, 3000])
ug = np.array([-10., -4.61])
vg = np.array([0, 0])

zsubs = np.array([0, 1500, 2100, 3000])
w = np.array([0, -0.65, 0., 0.])

zrad = np.array([0, 1500, 2500, 3000])
Qr = np.array([-2.0, -2.0, 0, 0])

zadv = np.array([0, 300, 500, 3000])
dqtdtadv = np.array([-1.2, -1.2, 0, 0])

# Plot the initial conditions
print 'Plotting the initial conditions'
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 2, 1)
ax1.plot(thetal, zthermo, 'k')
ax1.set_xlabel('$\\theta_{l}$ (K)')
ax1.set_ylabel('height (m)')
ax1.set_xticks([300, 304, 308])
ax1.set_xlim([298, 310])
ax1.set_ylim([0, 2500])
ax1.fill_between(x = [298, 310], y1 = [500, 500], y2 = [1500, 1500], color = 'black', alpha = 0.25, edgecolor = 'none')
ax1.text(x = 307, y = 2100, s = '$\\theta_{l}$')

ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks([4, 8, 12, 16])
ax2.set_xlim([2, 18])
ax2.plot(qt, zthermo, 'k')
ax2.set_xlabel('q$_{T}$ (g kg$^{-1}$)')
ax2.set_ylim([0, 2500])
ax2.text(x = 5, y = 2100, s = 'q$_{T}$')

ax3 = fig1.add_subplot(1, 2, 2)
ax3.plot(u, zwind, 'k')
ax3.set_xlim([-10, 2])
ax3.set_ylim([0, 2500])
ax3.set_xlabel('U (m s$^{-1}$)')
ax3.fill_between(x = [-10, 2], y1 = [500, 500], y2 = [1500, 1500], color = 'black', alpha = 0.25, edgecolor = 'none')
ax3.set_xticks([-8, -4, 0])
ax3.text(x = -8, y = 2100, s = 'U')

ax4 = ax3.twiny()
ax4.plot(v, zwind, 'k--')
ax4.set_xlim([-10, 2])
ax4.set_ylim([0, 2500])
ax4.set_xlabel('V (m s$^{-1}$)')
ax4.set_xticks([-8, -4, 0])
ax4.text(x = -1, y = 2100, s = 'V')
plt.savefig('../BOMEX_figs/initial_conditions.png', dpi = 150)
plt.close('all')

# Plot the large-scale forcings
print 'Plotting the large-scale forcings'
fig2 = plt.figure(tight_layout = True)
ax1 = fig2.add_subplot(1, 3, 1)
ax1.plot(ug, zgeo, 'k')
ax1.set_xlim([-10, 2])
ax1.set_xticks([-8, -4, 0])
ax1.set_ylim([0, 2500])
ax1.set_xlabel('U$_{g}$ (m s${-1}$)')
ax1.set_ylabel('height (m)')
ax1.text(x = -8, y = 2100, s = 'U$_{g}$')
ax1.fill_between(x = [-10, 2], y1 = [500, 500], y2 = [1500, 1500], color = 'black', alpha = 0.25, edgecolor = 'none')

ax2 = ax1.twiny()
ax2.plot(vg, zgeo, 'k--')
ax2.set_xlim([-10, 2])
ax2.set_xticks([-8, -4, 0])
ax2.set_ylim([0, 2500])
ax2.set_xlabel('V$_{g}$ (m s$^{-1}$)')
ax2.text(x = -2, y = 2100, s = 'V$_{g}$')

ax3 = fig2.add_subplot(1, 3, 2)
ax3.plot(w, zsubs, 'k')
ax3.set_ylim([0, 2500])
ax3.set_xlim([-0.8, 0.2])
ax3.set_xticks([-0.6, -0.3, 0])
ax3.set_xlabel('W (cm s$^{-1}$)')
ax3.fill_between(x = [-0.8, 0.2], y1 = [500, 500], y2 = [1500, 1500], color = 'black', alpha = 0.25, edgecolor = 'none')

ax4 = fig2.add_subplot(1, 3, 3)
ax4.plot(Qr, zrad, 'k')
ax4.set_xlim([-2.5, 0.5])
ax4.set_xticks([-2, -1, 0])
ax4.set_ylim([0, 2500])
ax4.set_xlabel('Q$_{r}$ (K day$^{-1}$)')
ax4.text(x = -1.5, y = 2100, s = 'Q$_{r}$')
ax4.fill_between(x = [-2.5, 0.5], y1 = [500, 500], y2 = [1500, 1500], color = 'black', alpha = 0.25, edgecolor = 'none')

ax5 = ax4.twiny()
ax5.plot(dqtdtadv, zadv, 'k--')
ax5.set_xlim([-2.5, 0.5])
ax5.set_xticks([-2., -1., 0])
ax5.set_ylim([0, 2500])
ax5.set_xlabel('$\\frac{\partial q_{T}}{\partial t}_{adv}$ (10$^{-8}$ kg kg$^{-1}$ s$^{-1}$)')
ax5.text(x = -0.5, y = 2100, s = '$\\frac{\partial q_{T}}{\partial t}$')

plt.savefig('../BOMEX_figs/largescale_forcings.png', dpi = 150)
plt.close('all')


