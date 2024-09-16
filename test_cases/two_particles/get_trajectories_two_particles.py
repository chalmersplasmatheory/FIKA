# Import libraries
import numpy as np
from numpy import sin, sqrt, cos
from scipy.integrate import cumulative_trapezoid
from scipy.constants import c, m_e, epsilon_0, e, hbar
import h5py

# Indicate parameters for the trajectory calculation
dt_ret   = 0.01e-15                                       # timestep
Tsim     = 4e-12                                          # total simulation time
E_MeV    = 150.0                                          # electron energy
r_beta   = 1.2e-6                                         # oscillation amplitude
gamma    = E_MeV / 0.51099895                             # gamma factor
t_ret    = np.arange(0, Tsim, dt_ret)                     # time array
ne       = 5e24                                           # plasma electron background density
omega_p  = sqrt(e * e * ne/(m_e * epsilon_0))             # plasma electron frequency
omega_b  = omega_p/sqrt(2 * gamma)                        # betatron frequency

# Calculate trajectory and corresponding velocities
y        = -r_beta * sin(omega_b * t_ret)                 # y positions electrons
vy       = np.zeros(np.size(y))                           # prealocate y velocity of electrons

for i in range(1,len(vy)):
  vy[i]  = (y[i]-y[i-1])/dt_ret                           # calculate velocity from y positions

vy[0]    = vy[1]                                              # change the initial velocity to next time step velocity, avoid 0 value  
by       = vy / c                                             # normalized y velocity
vec_ones = np.ones(np.size(by))                               # unit vector in the size of by
bx       = np.sqrt(vec_ones-vec_ones/(gamma*gamma)-by*by)     # normalized x velocity from gamma and by
vx       = bx * c                                             # x velocity
x        = cumulative_trapezoid(vx, t_ret, initial=t_ret[0])  # get x positions from time - integration of vx
z        = np.zeros(np.size(x))                               # z positions - empty array (no motion)
vz       = np.zeros(np.size(x))                               # z velocity - empty array (no motion)

# Prepare dataset for 1 particle for HDF output
ids = [1]
data = {
    1: {
        'time': t_ret,
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,},
}

# Write the dataset for one particle into a .h5 file
with h5py.File('test_cases/two_particles/one/particle_trajectory.h5', 'w') as f:
    for id, arrays in data.items():
        grp = f.create_group(str(id))
        for key, values in arrays.items():
            grp.create_dataset(key, data = values)

# Print update after the file creation
print("HDF5 file 'test_cases/two_particles/one/particle_trajectory.h5' created successfully.")

# Prepare datasets for 2 identical particles for HDF output
ids = [1,2]
data = {
    1: {
        'time': t_ret,
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,},
    2: {
        'time': t_ret,
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,
    },
}

# Write the dataset for two particles into a .h5 file
with h5py.File('test_cases/two_particles/two/particle_trajectories.h5', 'w') as f:
    for id, arrays in data.items():
        grp = f.create_group(str(id))
        for key, values in arrays.items():
            grp.create_dataset(key, data = values)

# Print update after the file creation
print("HDF5 file 'test_cases/two_particles/two/particle_trajectories.h5' created successfully.")
