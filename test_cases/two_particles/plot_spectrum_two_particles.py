# Import matplotlib for plotting and simulation parameters from get_trajectories_two_particles
import matplotlib
import matplotlib.pyplot as plt
from get_trajectories_two_particles import *

# Folder fit the FIKA output
folder     = 'test_cases/two_particles/'

# Read spectrum for one particle
with h5py.File(folder + 'one/final_spectrum.h5', 'r') as file:
    ene_1        = np.array(file['ene'])
    spectrum_1   = np.array(file['spectrum_ene'])
    t_1          = np.array(file['t'])
    spectrum_t_1 = np.array(file['spectrum_t'])

# Read spectrum for two indentical particles
with h5py.File(folder + 'two/final_spectrum.h5', 'r') as file:
    ene_2        = np.array(file['ene'])
    spectrum_2   = np.array(file['spectrum_ene'])
    t_2          = np.array(file['t'])
    spectrum_t_2 = np.array(file['spectrum_t'])

# Initiate the plot of spectrum at the observer position. 
# Parameters r = 1, phi = 0, theta = pi/2 were used in input.py.
# Other parameters used were E_slice_eV = 100, t_slice_s = 1e-19, E_radMax_eV = 25000.
fig, ax = plt.subplots()
# Plot spectrum of two same electrons, energies are in keV
plt.plot(ene_2/1000, spectrum_2 / 1e-15, label = r'2 $e^-$',color='black') 
# Put ticks inside the box on every axis
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
# Put axis labels
plt.xlabel('$E_{ph}$ [keV]', fontsize=14)
plt.ylabel(r'$\frac{d^2 W}{d\omega d\Omega}$ [$10^{-15}$ $\rm{J^{-1}}$ eV $\rm{sr}^-1$]', fontsize=14)
# Set the tick fontsize
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
# Set x limit for visualization
plt.xlim([0,25])
# Plot spectrum of one electron
plt.plot(ene_1/ 1000, spectrum_1 / 1e-15, label = '1 $e^-$', color='orange') 
# Set the plot title
plt.title(r' $E_e = $ 150 MeV, $r_{\beta}$ = 1.2 $\rm{\mu m}$',fontsize=14)
# Put the legend to the figure
plt.legend(frameon = False, fontsize = 14)
# Fit the figure to the box
plt.tight_layout()
# Save the figure
plt.savefig('test_cases/two_particles/energy_spectrum.png')

# Prepare similar figure for the temporal profile and save it
fig, ax = plt.subplots()
plt.plot(t_2/1e-15, spectrum_t_2, label = r'2 $e^-$',color='black') 
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
plt.xlabel('$E_{ph}$ [keV]', fontsize=14)
plt.ylabel(r'$\frac{d^2 W}{d t d\Omega}$ [$\rm{J^{-1}}$ $\rm{s^{-1}}$ $\rm{sr}^-1$]', fontsize=14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(t_1/ 1e-15, spectrum_t_1 , label = '1 $e^-$', color='orange') 
plt.title(r' $E_e = $ 150 MeV, $r_{\beta}$ = 1.2 $\rm{\mu m}$',fontsize=14)
plt.legend(frameon = False, fontsize = 14)
mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
plt.gca().yaxis.set_major_formatter(mf)
ax.yaxis.offsetText.set_fontsize(14)
plt.tight_layout()
plt.savefig('test_cases/two_particles/temporal_profile.png')

# Show plots
plt.show()