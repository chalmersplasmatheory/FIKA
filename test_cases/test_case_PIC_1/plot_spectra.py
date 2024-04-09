# This is an example of electron trajectories obtained from a particle-in-cell simulation with the SMILEI code.
# The radiation has synchrotron-like character, corresponding to betatron radiation from laser wakefield acceleration.
# Electron beam is moving in the positive x-direction undergoing small transverse betatron oscillations in the y & z directions.
# Here, the macroparticle spectra (PIC_macroparticle_weights = False) and real-particle spectrum (PIC_macroparticle_weights = True) are plotted.
# The options in the input file were as follows:
# The position of the observer was chosen 1 m distance from the beam on the x-axis (in the positive direction).
# It corresponds to r = 1, phi = 0, theta = pi/2.
# E_slice_eV = 10 and t_slice_s = 1e-16 were chosen in the spectral summation.
# For plotting, use python plot_spectra.py

# Import necessary libraries for data handling and plotting
import h5py
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# Function to plot the data
def plot_spectrum(x, y, x_label, y_label, title, color):
    fig, ax = plt.subplots()
    plt.plot(x, y, color=color)
    ax.tick_params(axis="x", direction="in", bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis="y", direction="in", bottom=True, top=True, left=True, right=True)
    plt.title(title, fontsize=14, color=color)
    plt.xlabel(x_label, fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
    plt.gca().yaxis.set_major_formatter(mf)
    ax.yaxis.offsetText.set_fontsize(14)

# Reading data from HDF5 file for the case without macroparticle weights
folder = 'test_cases/test_case_PIC_1/without_macroparticle_weights/'
with h5py.File(folder + 'final_spectrum.h5', 'r') as file:
  ene          = np.array(file['ene'])
  spectrum_ene = np.array(file['spectrum_ene'])
  t            = np.array(file['t'])
  spectrum_t   = np.array(file['spectrum_t'])

# Reading data from HDF5 file for the case with macroparticle weights
folder_w = 'test_cases/test_case_PIC_1/with_macroparticle_weights/'
with h5py.File(folder_w + 'final_spectrum.h5', 'r') as file_w:
    ene_w          = np.array(file_w['ene'])
    spectrum_ene_w = np.array(file_w['spectrum_ene'])
    t_w            = np.array(file_w['t'])
    spectrum_t_w   = np.array(file_w['spectrum_t'])

# Plotting the energy spectrum without macroparticle weights, energy is in keV
plot_spectrum(ene / 1000, spectrum_ene,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{dE d\Omega}$ [$\rm{J}$ $\rm{eV^{-1}}$ $\rm{sr}^-1$]', 'Without macroparticle weights', 'darkred')
plt.xlim([0,1])
plt.tight_layout()
plt.savefig('test_cases/test_case_PIC_1/macroparticles_ene_spectrum.png')

# Plotting the energy spectrum with macroparticle weights, energy is in keV
plot_spectrum(ene_w / 1000, spectrum_ene_w,'$E_{ph}$ [keV]',r'$\frac{d^2 W}{dE d\Omega}$ [$\rm{J}$ $\rm{eV^{-1}}$ $\rm{sr}^-1$]', 'With macroparticle weights', 'green')
plt.xlim([0,1])
plt.tight_layout()
plt.savefig('test_cases/test_case_PIC_1/real_particles_ene_spectrum.png')

# Plotting the temporal profile without macroparticle weights, seconds are converted to fs
plot_spectrum(t / 1e-15, spectrum_t * 1e-15,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{fs}^{-1}$ $\rm{sr}^-1$]', 'Without macroparticle weights', 'darkred')
plt.savefig('test_cases/test_case_PIC_1/macroparticles_spectrum_t.png')
plt.tight_layout()

# Plotting the temporal profile with macroparticle weights, seconds are converted to fs
plot_spectrum(t_w / 1e-15, spectrum_t_w * 1e-15,'$t$ [fs]', r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{fs}^{-1}$ $\rm{sr}^-1$]', 'With macroparticle weights', 'green')
plt.tight_layout()
plt.savefig('test_cases/test_case_PIC_1/real_particles_spectrum_t.png')

plt.show()

