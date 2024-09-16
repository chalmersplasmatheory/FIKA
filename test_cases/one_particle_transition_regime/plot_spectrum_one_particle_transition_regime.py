# Import matplotlib library for plotting and simulation parameters
import matplotlib.pyplot as plt
import matplotlib
from get_trajectory_one_particle_transition_regime import *

# Read output of the FIKA simulation
with h5py.File('test_cases/one_particle_transition_regime/final_spectrum.h5', 'r') as file:
  ene          = np.array(file['ene'])
  spectrum_ene = np.array(file['spectrum_ene'])
  signal_time  = np.array(file['t'])
  spectrum_t   = np.array(file['spectrum_t'])

# Plot temporal profile at the observer position. Parameters r = 1, phi = 0, theta = pi/2 were used in input.py.
# Other parameters used were E_slice_eV = 0.1, t_slice_s = 1e-18.
fig, ax = plt.subplots()
plt.title(r'$E_{e} = $ 25 MeV, $r_{\beta}$ = 0.5 $\rm{\mu m}$',fontsize=14)
plt.plot(signal_time / 1e-15, spectrum_t)
plt.xlabel('$t$ [fs]', fontsize = 14)
plt.ylabel(r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]', fontsize=14)
mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
plt.gca().yaxis.set_major_formatter(mf)
ax.yaxis.offsetText.set_fontsize(14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
plt.tight_layout()
plt.savefig('test_cases/one_particle_transition_regime/temporal_profile.png')

# Plot energy spectrum
fig, ax = plt.subplots()
plt.plot(ene, spectrum_ene) 
plt.title(r'$E_{e} = $ 25 MeV, $r_{\beta}$ = 0.5 $\rm{\mu m}$',fontsize=14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True) 
plt.xlabel('$E_{ph}$ [eV]', fontsize=14)
plt.ylabel(r'$\frac{d^2 W}{d\omega d\Omega}$ [ $\rm{J}$ $\rm{s}$ $\rm{rad}^{-1}$ $\rm{sr}^{-1}$]', fontsize=14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim([0,300])
plt.tight_layout()
plt.savefig('test_cases/one_particle_transition_regime/energy_spectrum.png')

# Calculate expected fundamental frequency with analytical formula and compare it with the code outcome.
K       = 1.33e-10 * sqrt(gamma*n_e/1e6) * r_beta/1e-6
omega_1 = 2 * gamma**2 * omega_beta/ (1 + K**2 / 2)
max_s   = np.argmax(spectrum_ene)
print('Expected fundamental frequency according to formula 2*gamma^2*omega_beta/(1+K^2/2), where K = 1.33*10^-10*r_beta*(gamma*n_e)^0.5 is', int(round(omega_1*hbar/e)), 'eV.')
print('The fundamental frequency in the spectrum calculated by FIKA is ' +str(int(round(ene[max_s])))+' eV.')
plt.show()
