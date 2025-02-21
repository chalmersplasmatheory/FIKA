# Import libraries and data info from get_trajectory_one_particle_undulator script
import matplotlib.pyplot as plt
import matplotlib
from scipy.special import kv
from get_trajectory_circle import *
import os
import sys

# Create a string representing the closest common directory
base_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(base_folder, "../../"))
from input import *

# Read data from FIKA simulation:
                #/Users/alberthansson/Documents/Shalmerc/FIKA/GITHUB/INTENSE/test_cases/circle/final_spectrum.h5
with h5py.File(base_folder + '/test_cases/one_particle_circular_orbit/final_spectrum.h5', 'r') as file:
    ene          = np.array(file['ene'])
    spectrum_ene = np.array(file['spectrum_ene'])
    signal_time  = np.array(file['t'])
    spectrum_t   = np.array(file['spectrum_t'])

# Plot temporal profile at the observer position. Parameters r = 1, phi = 0, theta = pi/2 were used in input.py.
fig, ax = plt.subplots()
# Write important info into the title
plt.title(r'$E_{e} = $10 MeV, $r_o$ = 0.05 $\rm{\mu m}$', fontsize=14)
# Plot temporal profile, time is in the units of fs
plt.plot(signal_time / 1e-15, spectrum_t)
# Put x & y labels
plt.xlabel('$t$ [fs]', fontsize = 14)
plt.ylabel(r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]', fontsize=14)
#Format ticks fontsize, position and orientation
mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
plt.gca().yaxis.set_major_formatter(mf)
ax.yaxis.offsetText.set_fontsize(14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
# Fit figure to the figure box and save
plt.tight_layout()
plt.savefig(base_folder + '/test_cases/one_particle_circular_orbit/temporal_profile.png')
plt.savefig(base_folder + '/test_cases/one_particle_circular_orbit/temporal_profile.pdf')
plt.show()

# Plot energy spectrum
fig, ax = plt.subplots()

plt.plot(ene/1000, spectrum_ene*hbar / max(spectrum_ene*hbar), alpha = 0.3, color = 'blue')
plt.title(r'$E_{e} = $10 MeV, $r_o$ = 0.05 $\rm{\mu m}$',fontsize=14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True) 
plt.xlabel('$E_{ph}$ [keV]', fontsize=14)
plt.ylabel(r'$\frac{d^2 W}{d\omega d\Omega} / \left(\frac{d^2 W}{d\omega d\Omega}\right)_{max}$', fontsize=14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim(0,3e2)
plt.tight_layout()


theta = pi/2 - theta
omega = np.linspace(0, (ene[-1]*e/hbar), 10000)
omega_c = 3/2 * gamma**3 * c / r_beta
xi = omega / omega_c * (1 + (gamma*theta)**2)**1.5/2
test = e**2*gamma*2/(16*pi**3*epsilon_0*c) * (omega/omega_c)**2 * (1 + (gamma*theta)**2)**2
dWdwdo = e**2*gamma*2/(16*pi**3*epsilon_0*c) * (omega/omega_c)**2 * (1 + (gamma*theta)**2)**2 * ( kv(2/3, xi)**2 + (gamma * theta)**2 / (1 + (gamma*theta)**2) * kv(1/3, xi)**2 )


plt.plot(omega*hbar/e / 1000, dWdwdo/hbar / max(dWdwdo[1:]/hbar), color = 'orange')


energy_guess_interval = [1, len(omega) - 1]
while np.abs(energy_guess_interval[0] - energy_guess_interval[1]) > 1:
    new_guess = int((energy_guess_interval[0] + energy_guess_interval[1])/2)
    If = np.sum(dWdwdo[1:new_guess])
    Ie = np.sum(dWdwdo[new_guess:-1])
    if (If > Ie):
        energy_guess_interval = [energy_guess_interval[0], new_guess]
    else:
        energy_guess_interval = [new_guess, energy_guess_interval[1]]
critical_energy_theory = omega[energy_guess_interval[0]]*hbar/e # eV
vertical_limit = dWdwdo[energy_guess_interval[0]] / max(dWdwdo[1:])

energy_guess_interval = [0, len(ene) - 1]
while np.abs(energy_guess_interval[0] - energy_guess_interval[1]) > 1:
    new_guess = int((energy_guess_interval[0] + energy_guess_interval[1])/2)
    If = np.sum(spectrum_ene[0:new_guess])
    Ie = np.sum(spectrum_ene[new_guess:-1])
    if (If > Ie):
        energy_guess_interval = [energy_guess_interval[0], new_guess]
    else:
        energy_guess_interval = [new_guess, energy_guess_interval[1]]
critical_energy = ene[energy_guess_interval[0]] # eV
plt.plot(np.array([critical_energy, critical_energy]) / 1000,[0, vertical_limit], zorder = 10, color = 'blue', linestyle = '--')
plt.plot(np.array([critical_energy_theory, critical_energy_theory]) / 1000,[0, vertical_limit], zorder = 11, color = 'orange', linestyle = '--')



plt.xlim(0,3e2)
plt.legend(['FIKA calculation', 'Theoretical spectrum', r'$E_c$ calculated', r'$E_c$ theory'])
plt.grid()
plt.savefig(base_folder + '/test_cases/one_particle_circular_orbit/energy_spectrum.png')
plt.savefig(base_folder + '/test_cases/one_particle_circular_orbit/energy_spectrum.pdf')
plt.show()
