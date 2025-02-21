# Import libraries
import numpy as np
import os
from numpy import cos,sin,sqrt
from matplotlib import pylab as plt
from scipy.constants import c, m_e, epsilon_0, e, hbar, pi
import h5py

# Create a string representing the closest common directory
base_folder = os.path.dirname(os.path.abspath(__file__))

# Indicate parameters for the trajectory calculation
dt_ret     = 0.01e-16                            # timestep
E_MeV      = 10                                # total electron energy (kinetic + rest) in MeV
r_beta     = 0.05e-6                             # orbit radius
gamma      = E_MeV / 0.51099895                  # gamma factor
omega_beta = c/r_beta * sqrt(1-1/(gamma**2))     # angular frequency of orbit
t_1_lap = 2*pi/omega_beta                        # duration of one lap
Tsim       = 5*t_1_lap                          # total simulation time
t_ret      = np.arange(0, Tsim + 0, dt_ret)      # time array

# Calculate trajectory and corresponding velocities
y          = r_beta * sin(omega_beta * t_ret)              # y positions of the electron
vy         = omega_beta * r_beta * cos(omega_beta * t_ret) # y velocity of electron
ay = (-1) * omega_beta**2 * r_beta * sin(omega_beta * t_ret)
by       = vy / c                                          # normalized y velocity
vec_ones = np.ones(np.size(by))                            # unit vector in the size of by
bx       = vy / c                                          # normalized x velocity calculated from gamma and by

x        = r_beta * cos(omega_beta * t_ret)                # x positions of the electron
vx       = (-1) * omega_beta * r_beta * sin(omega_beta * t_ret)   # x velocity
ax = (-1) * omega_beta**2 * r_beta * cos(omega_beta * t_ret)
z        = np.zeros(np.size(x))                            # z positions - zero array (no motion in z)
vz       = np.zeros(np.size(x))                            # z velocity - zero array (no motion in z)

a = omega_beta * omega_beta * r_beta / c
beta = sqrt(vx**2 + vy**2)[0]/c
#print(sqrt(vx**2 + vy**2)[0]/c)
# gamma = 1/sqrt(1-beta**2)
#print(beta, gamma)
w_ph = 3/2 * gamma**3 * c / r_beta
E_ph = w_ph * hbar / e
#print(w_ph, E_ph)


plt.figure(figsize=(4,4))
plt.title('In-plane view of orbital motion')
plt.plot(x / 1e-9, y / 1e-9)
plt.xlim(-2.1*r_beta / 1e-9, 2.1*r_beta / 1e-9)
plt.ylim(-2.1*r_beta / 1e-9, 2.1*r_beta / 1e-9)
plt.xlabel('x position [nm]')
plt.ylabel('y position [nm]')
plt.show()

plt.figure(figsize=(4,4))
plt.title('Temporal profile of orbital ocsillation')
plt.plot(t_ret / 1e-15, y / 1e-9)
plt.xlabel('time [fs]')
plt.ylabel('y position [nm]')
plt.show()


# Prepare dataset for HDF output
ids = [1]
data = {
    1: {
        'time': t_ret,
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,
    },
}

# Write the dataset for one particle into a .h5 file

with h5py.File(base_folder + '/particle_trajectory.h5', 'w') as f:
    for id, arrays in data.items():
        grp = f.create_group(str(id))
        for key, values in arrays.items():
            grp.create_dataset(key, data=values)

# Print status when all done
print("HDF5 file 'particle_trajectory.h5' created successfully.")
