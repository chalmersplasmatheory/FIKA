import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from input import *
import os
from utils import cartesian_interpolate

spatial_profile_path = os.path.join(output_folder, 'spatial_profile.txt')

# Read data from saved file
with open(spatial_profile_path, 'r') as file:
    all_data = np.array([float(line.strip()) for line in file.readlines()[1:]])
    plot_settings = all_data[:3]
    data_list = all_data[3:]
data_grid = data_list.reshape((int(np.sqrt(len(data_list))), int(np.sqrt(len(data_list)))))

# Set up the deflection grid
resolution_s = int(plot_settings[0])
boundaries_s = [float(plot_settings[1]), float(plot_settings[2])]
phi_angles = np.linspace(boundaries_s[0], boundaries_s[1], resolution_s)
theta_angles = phi_angles.copy()
phi_grid, theta_grid = np.meshgrid(phi_angles, theta_angles, indexing='ij')
angles = np.stack([phi_grid.ravel(), theta_grid.ravel()], axis=-1)

# Compute total energy
dphi = (boundaries_s[1] - boundaries_s[0]) / (resolution_s - 1)
dtheta = dphi
total_diff_energy = np.sum(data_grid)
print(f"Total energy: {total_diff_energy * dphi * dtheta} J")

# Plot
extent = [boundaries_s[0]*1000, boundaries_s[1]*1000, boundaries_s[0]*1000, boundaries_s[1]*1000]
fig = plt.figure(figsize=(7, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.05], wspace=0.05)
ax = fig.add_subplot(gs[0])
PHI, THETA = np.meshgrid(phi_angles, phi_angles)

if do_interpolation:
    PHI, THETA, data_grid = cartesian_interpolate(data_grid, PHI, THETA)

im = ax.pcolormesh(PHI * 1000, THETA * 1000, data_grid.T)

ax.set_aspect('equal', adjustable='box')
cax = fig.add_subplot(gs[1])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label(r"$\text{d}W/\text{d}\Omega$", fontsize=14)
ax.set_xlabel(r'$\varphi$ [mrad]', fontsize=14)
ax.set_ylabel(r'$\theta$ [mrad]', fontsize=14)
ax.set_title("Optimal case", fontsize=16)
ax.tick_params(labelsize=12)
ticks = np.round(np.linspace(np.min(phi_angles), np.max(phi_angles), 7), 3) * 1000
ax.set_xticks(ticks)
ax.set_yticks(ticks)

plt.tight_layout()
plt.show()