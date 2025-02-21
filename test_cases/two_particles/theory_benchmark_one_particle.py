# Import matplotlib for plotting and simulation parameters from get_trajectories_two_particles
import matplotlib
import matplotlib.pyplot as plt
from get_trajectories_two_particles import *
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import pandas as pd
def round_to_odd(f):
  return int(np.ceil(f) // 2 * 2 + 1)
 
 import os

# Create a string representing the closest common directory
base_folder = os.path.dirname(os.path.abspath(__file__))

# Read spectrum for one particle
with h5py.File(base_folder + '/one/final_spectrum.h5', 'r') as file:
    ene_1        = np.array(file['ene'])/hbar*e
    spectrum_1   = np.array(file['spectrum_ene'])/(hbar/e)
    t_1          = np.array(file['t'])
    spectrum_t_1 = np.array(file['spectrum_t'])

# Read spectrum for two indentical particles
with h5py.File(base_folder + '/two/final_spectrum.h5', 'r') as file:
    ene_2        = np.array(file['ene'])
    spectrum_2   = np.array(file['spectrum_ene'])
    t_2          = np.array(file['t'])
    spectrum_t_2 = np.array(file['spectrum_t'])

# Initiate the plot of spectrum at the observer position. 
# Parameters r = 1, phi = 0, theta = pi/2 were used in input.py.
# Other parameters used were E_slice_eV = 100, t_slice_s = 1e-19.
fig, ax = plt.subplots()

# Plot spectrum of one electron
#plt.plot(ene_1, spectrum_1, label = '1 $e^-$', color='orange')
#plt.xlim(0,1e20) 
# Set the plot title
plt.title(r' $E_e = $ 150 MeV, $r_{\beta}$ = 1.2 $\rm{\mu m}$',fontsize=14)
# Put the legend to the figure
plt.legend(frameon = False, fontsize = 14)
# Fit the figure to the box
plt.tight_layout()
# Save the figure
plt.savefig(base_folder + '/energy_spectrum.png')


# Window size
window_size = 500


# Create a window for convolution that represents the moving average
window = np.ones(window_size) / window_size

# Calculate the moving average using convolution
moving_average = np.convolve(spectrum_1, window, 'same')



cumulative_spectrum = np.cumsum(moving_average)
# Finding the point corresponding to the half the power radiated
# The spectrum on the axis is *1.54 higher than the whole beam critical energy
halfway_point       = cumulative_spectrum[-1] / 2 /1.54  
index_halfway       = np.searchsorted(cumulative_spectrum, halfway_point)
energy_crit_code    = ene_1[index_halfway] 
print('Critical energy of the calculated spectrum of one particle is '+ str(energy_crit_code)+'.')

#plt.plot(ene_1/1000,moving_average)
index_of_max = np.argmax(moving_average)

print("Index of maximum value:", ene_1[index_of_max])

k_p = omega_p/c
print(ene_1[index_halfway-1]/1000,ene_1[index_halfway+1]/1000,ene_1[index_halfway]/1000)
K                  = r_beta * k_p * sqrt(gamma/2.)
omega_crit_theory  = (3/2) * K * omega_b * gamma**2
energy_crit_theory =  omega_crit_theory  
print('Critical energy of one particle according to the theory is '+str(energy_crit_theory)+'.')
k_b = omega_b/c
a_b = k_b*r_beta*gamma
omm=3*a_b*gamma*gamma*omega_b
print(omega_crit_theory,omm, 5.21e-21*gamma*gamma*5e18*0.5/hbar*e)
print(ne/1e6,gamma, r_beta/1e-6)
nb=9
xii=np.linspace(0,10,100) #ene_1[1:100]/hbar*e/omega_crit_theory
from scipy.special import kv
dii=xii*xii*(kv(2/3, xii)**2)#*nb*6*e*e/(c*np.pi*np.pi)*gamma*gamma
#plt.plot(xii,dii,color='green')
print(np.nanmax(dii))
#plt.xlim([0,5e19])
#NORMALIZED
plt.plot(ene_1/omm, spectrum_1/max(spectrum_1), color='grey',label='code')
plt.plot(ene_1/omm, moving_average/max(moving_average), color='black',label =' code moving average')
plt.plot(xii,dii/np.nanmax(dii),color='red',label='theory')
plt.legend()
plt.xlim([0,6])
'''
# Prepare similar figure for the temporal profile and save it
fig, ax = plt.subplots()
#plt.plot(t_2/1e-15, spectrum_t_2, label = r'2 $e^-$',color='black') 
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
plt.xlabel('$t$ [fs]', fontsize=14)
plt.ylabel(r'$\frac{d^2 W}{d t d\Omega}$ [$\rm{J^{-1}}$ $\rm{s^{-1}}$ $\rm{sr}^{-1}$]', fontsize=14)
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
'''
# Show plots
plt.show()
