import matplotlib.pyplot as plt
import matplotlib
from get_trajectories_one_particle_undulator import *
from scipy.signal import savgol_filter
with h5py.File('test_cases/one_particle_undulator/final_spectrum.h5', 'r') as file:
    ene        = np.array(file['ene'])
    spectrum    = np.array(file['spectrum_ene'])
    signal_time = np.array(file['t'])
    signal_y    = np.array(file['spectrum_t'])

def round_to_odd(f):
  return int(np.ceil(f) // 2 * 2 + 1)

fig, ax = plt.subplots()
plt.title(r'$E_{tot}} = $ 10 MeV, $r_{\beta}$ = 0.05 $\rm{\mu m}$',fontsize=14)
plt.plot(signal_time / 1e-15, signal_y)
plt.xlabel('$t$ [fs]', fontsize = 14)
plt.ylabel(r'$\frac{d^2 W}{\rm{d}t \rm{d}\Omega}$ [$\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^-1$]', fontsize=14)

#spectrum_smoothed_w = savgol_filter(signal_y, round_to_odd((len(signal_y))/50), 3)
#plt.plot(signal_time / 1e-15, spectrum_smoothed_w, color = 'darkblue',linestyle = 'dashed', label = 'smoothed data weighted') 
mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
plt.gca().yaxis.set_major_formatter(mf)
ax.yaxis.offsetText.set_fontsize(14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True)
plt.tight_layout()

#plt.savefig('figs/test_case_1_temporal_profile.pdf')

K       = 1.33e-10 * sqrt(gamma*n_e) * r_beta
omega_1 = 2 * gamma**2 * omega_beta/ (1 + K**2 / 2)
print('Expected fundamental frequency according to formula 2*gamma^2*omega_beta/1+K^2/2, where K = 1.33*10^-10*r_beta*(gamma*n_e)^0.5 is ' + str(round(omega_1*hbar/e,1))+ 'eV')
max_s = np.argmax(spectrum)
print(round(ene[max_s],2))


fig, ax = plt.subplots()
plt.plot(ene, spectrum / 1e-21, label = 'raw data') 
plt.title(r'$E_{tot} = $ 10 MeV, $r_{\beta}$ = 0.05 $\rm{\mu m}$',fontsize=14)
ax.tick_params(axis="x",direction="in", bottom=True, top=True, left=True, right=True)
ax.tick_params(axis="y",direction="in", bottom=True, top=True, left=True, right=True) 
plt.xlabel('$E_{ph}$ [eV]', fontsize=14)
plt.ylabel(r'$\frac{d^2 I}{d\omega d\Omega}$ [$10^{-21}$ $\rm{J}$ $\rm{eV}^{-1}$ $\rm{sr}^-1$]', fontsize=14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim([0,50])
plt.ylim([0,10])
plt.tight_layout()
#plt.savefig('figs/test_case_1_spectrum.pdf')
#plt.savefig('figs/test_case_1_spectrum.png')
plt.show()

