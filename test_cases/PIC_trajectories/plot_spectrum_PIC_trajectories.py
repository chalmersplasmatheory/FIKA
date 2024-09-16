# The options in the input file were as follows: 
# The position of the observer was chosen 1 m distance from the beam on the x-axis (in the positive direction).
# It corresponds to r = 1, phi = 0, theta = pi/2. 
# Also, E_slice_eV = 10 and t_slice_s = 1e-19 were chosen.

# Import necessary libraries for data handling and plotting
import h5py
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.constants import c, epsilon_0, hbar, e

# Function to plot the data
def plot_spectrum(x, y, x_label, y_label, label, color, linestyle):
    plt.plot(x, y, color = color, linestyle = linestyle, label = label,  linewidth = 2)
    ax.tick_params(axis="x", direction="in", bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis="y", direction="in", bottom=True, top=True, left=True, right=True)
    plt.xlabel(x_label, fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
    plt.gca().yaxis.set_major_formatter(mf)
    ax.yaxis.offsetText.set_fontsize(14)

# Reading data from HDF5 file
def read_data_to_plot(folder):
  with h5py.File(folder + 'final_spectrum.h5', 'r') as file:
    ene          = np.array(file['ene'])
    spectrum_ene = np.array(file['spectrum_ene'])
    t            = np.array(file['t'])
    spectrum_t   = np.array(file['spectrum_t'])
    return ene, spectrum_ene, t, spectrum_t

# Reading data from HDF5 file for the case without macroparticle weights, serial calculation
ene_ser, spectrum_ene_ser, t_ser, spectrum_t_ser = read_data_to_plot('test_cases/PIC_trajectories/no_weights_serial/')
# Reading data from HDF5 file for the case with macroparticle weights, serial calculation
ene_w_ser, spectrum_ene_w_ser, t_w_ser, spectrum_t_w_ser = read_data_to_plot('test_cases/PIC_trajectories/weights_serial/')
# Reading data from HDF5 file for the case without macroparticle weights, parallel calculation on 4 mpi processes
ene_par, spectrum_ene_par, t_par, spectrum_t_par = read_data_to_plot('test_cases/PIC_trajectories/no_weights_parallel/')
# Reading data from HDF5 file for the case with macroparticle weights, parallel calculation on 4 mpi processes
ene_w_par, spectrum_ene_w_par, t_w_par, spectrum_t_w_par = read_data_to_plot('test_cases/PIC_trajectories/weights_parallel/')


# Plotting the energy spectrum without macroparticle weights, energy is in keV
fig, ax = plt.subplots()
plot_spectrum(ene_ser / 1000, spectrum_ene_ser ,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{d\omega d\Omega}$ [ $\rm{J}$ $\rm{s}$ $\rm{rad}^{-1}$ $\rm{sr}^{-1}$]',  'No weights, serial', 'pink', 'solid')
plot_spectrum(ene_par / 1000, spectrum_ene_par ,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{d\omega d\Omega}$ [ $\rm{J}$ $\rm{s}$ $\rm{rad}^{-1}$ $\rm{sr}^{-1}$]', 'No weights, parallel 4 mpi processes', 'red','dashed')
plt.xlim([0,3])
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.savefig('test_cases/PIC_trajectories/macroparticles_ene_spectrum.png')

# Plotting the energy spectrum with macroparticle weights, energy is in keV
fig, ax = plt.subplots()
plot_spectrum(ene_w_ser / 1000, spectrum_ene_w_ser,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{d\omega d\Omega}$ [ $\rm{J}$ $\rm{s}$ $\rm{rad}^{-1}$ $\rm{sr}^{-1}$]', 'Weights, serial', 'palegreen', 'solid')
plot_spectrum(ene_w_par / 1000, spectrum_ene_w_par ,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{d\omega d\Omega}$ [ $\rm{J}$ $\rm{s}$ $\rm{rad}^{-1}$ $\rm{sr}^{-1}$]', 'Weights,  parallel 4 mpi processes', 'green', 'dashed')
plt.xlim([0,3])
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.savefig('test_cases/PIC_trajectories/real_particles_ene_spectrum.png')

# Plotting the temporal profile without macroparticle weights, seconds are converted to fs
fig, ax = plt.subplots()
plot_spectrum(t_ser / 1e-15, spectrum_t_ser ,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]',  'No weights, serial', 'pink', 'solid')
plot_spectrum(t_par / 1e-15, spectrum_t_par ,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]',  'No weights, parallel 4 mpi processes' ,'red','dashed')
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.savefig('test_cases/PIC_trajectories/macroparticles_spectrum_t.png')

# Plotting the temporal profile with macroparticle weights, seconds are converted to fs
fig, ax = plt.subplots()
plot_spectrum(t_w_ser / 1e-15, spectrum_t_w_ser ,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]',  'Weights, serial', 'palegreen', 'solid')
plot_spectrum(t_w_par / 1e-15, spectrum_t_w_par,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]',  'Weights,  parallel 4 mpi processes', 'green', 'dashed')
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.savefig('test_cases/PIC_trajectories/real_particles_spectrum_t.png')

plt.show()

