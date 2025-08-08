import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
from scipy.signal import savgol_filter
base_folder = os.path.dirname(os.path.abspath(__file__))

# Definitions of folders and configuration
output_file_name = 'multiplot' 
case_name = 'PIC_multicase' # name of the folder in which you keep the case folders
sim_folders     = ['case_1','case_2','case_3','case_4'] # names of the case folders
legend_names= ['L1', 'L2', 'L3', 'L4']
colors      = ['black', 'green', 'red', 'blue'] # this color list must hold at least as many colors as there are cases


file_types = [
    ('Spectral profile', r'$\frac{\text{d}^2 W}{\text{d}\omega \text{d}\Omega} |_{\text{on-axis}}$ [ $\rm{J}$ $\rm{s}$ $\rm{sr}^{-1}$]'),
    ('Temporal profile', r'$\frac{\text{d}^2 W}{\text{d}t \text{d}\Omega} |_{\text{on-axis}} $ [ $\rm{J}$ $\rm{s}^{-1}$ $\rm{sr}^{-1}$]')
]

# Store the handles and labels
handles_list = []
labels_list = []

# Preparation of data storage
number_of_simulations = len(sim_folders)
critical_energies = ['Critical energy [eV]']
pulse_times = ['Pulse time [as]']
I_peaks = ['Intesity peak [J/sr]']

fig_t, ax_t = plt.subplots(figsize=(8, 5))
fig_e, ax_e = plt.subplots(figsize=(8, 5))

# Loop through all cases
for i, folder in enumerate(sim_folders):
    file_path = base_folder + f'/test_cases/{case_name}/{folder}/final_spectrum.h5'

    with h5py.File(file_path, 'r') as file:
        ene          = np.array(file['ene'])
        spectrum_ene = np.array(file['spectrum_ene'])
        signal_time  = np.array(file['t'])
        spectrum_t   = np.array(file['spectrum_t'])

    color = colors[i % len(colors)]


    # Plotting the temporal and spectral profiles 
    savgol_spectrum_t = savgol_filter(spectrum_t[::1], int(len(spectrum_t)/1000), 4)
    line, = ax_t.plot(signal_time*1e15, savgol_spectrum_t, linewidth=1, color = color)
    ax_t.set_xlabel('t [fs]', fontsize = 14)
    ax_t.set_title(file_types[1][0], fontsize = 16)

    line, = ax_e.plot(ene / 1000, spectrum_ene, linewidth=1, color=color)
    ax_e.set_xlabel('$E_{ph}$ [keV]', fontsize = 14)
    ax_e.set_xlim(0,10)
    ax_e.set_title(file_types[0][0], fontsize = 16)

    labels_list.append(legend_names[i])
    handles_list.append(line)


    # Algorithm for computing the critical energy
    energy_guess_interval = [0, len(ene) - 1]
    while np.abs(energy_guess_interval[0] - energy_guess_interval[1]) > 1:
        new_guess = int((energy_guess_interval[0] + energy_guess_interval[1])/2)
        If = np.sum(spectrum_ene[0:new_guess])
        Ie = np.sum(spectrum_ene[new_guess:-1])
        if (If > Ie):
            energy_guess_interval = [energy_guess_interval[0], new_guess]
        else:
            energy_guess_interval = [new_guess, energy_guess_interval[1]]
    critical_energy = energy_guess_interval[0]
    critical_energies.append((ene[critical_energy]))
    
    # Algorithm for computing the pulse time at full width half maximum
    HM = np.max(savgol_spectrum_t)*0.5
    for a in range(1, len(signal_time)-1):
        if savgol_spectrum_t[a] < HM and savgol_spectrum_t[a+1] > HM:
            tau_start_index = a
            break
    for a in range(1, len(signal_time)-1):
        if savgol_spectrum_t[a] > HM and savgol_spectrum_t[a+1] < HM:
            tau_end_index = a
            break
    
    pulse_duration_hm  = signal_time[tau_end_index]*1e18 - signal_time[tau_start_index]*1e18
    pulse_times.append((pulse_duration_hm)) # attoseconds

    # Storing the intensity peak 
    I_peak = np.max(savgol_spectrum_t)
    I_peaks.append((I_peak))


# Setting y axis parameters
ax_e.set_ylabel(file_types[0][1], fontsize=14)
ax_e.grid(True, alpha=0.35)
ax_e.tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=14)
ax_t.set_ylabel(file_types[1][1], fontsize=14)
ax_t.grid(True, alpha=0.35)
ax_t.tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=14)
fig_e.subplots_adjust(top=0.82)
fig_t.subplots_adjust(top=0.82)

# Adjust the spacing between subplots to ensure there is enough space
fig_t.legend(
    handles_list, 
    labels_list, 
    fontsize=14, 
    loc='upper center', 
    ncol=6,  # Set ncol to the max number of items expected in the longest row
    frameon=False, 
    handletextpad=0.5,  # Space between marker and label
    bbox_to_anchor=(0.53, 0.98),  # Adjust this slightly above the plot for centering
    columnspacing=1.5,  # Adjust columnspacing to control the appearance of spacing
    handlelength=1.5  # Adjust handle length for better spacing if needed
)
fig_e.legend(
    handles_list, 
    labels_list, 
    fontsize=14, 
    loc='upper center', 
    ncol=6,  # Set ncol to the max number of items expected in the longest row
    frameon=False, 
    handletextpad=0.5,  # Space between marker and label
    bbox_to_anchor=(0.53, 0.98),  # Adjust this slightly above the plot for centering
    columnspacing=1.5,  # Adjust columnspacing to control the appearance of spacing
    handlelength=1.5  # Adjust handle length for better spacing if needed
)

# Save and show the figure
plt.savefig(base_folder + '/test_cases/' + case_name + '/' + output_file_name + '.png', dpi=300)
plt.savefig(base_folder + '/test_cases/' + case_name + '/' + output_file_name + '.pdf', dpi=300)
plt.show()

# Produce text file with pulse characteristics
output_array = np.array([critical_energies, I_peaks, pulse_times]).T
np.savetxt(base_folder + '/test_cases/' + case_name + '/' + output_file_name + '_data.csv', output_array, fmt='%s', delimiter=',')
plt.rcParams['font.family'] = 'Helvetica'