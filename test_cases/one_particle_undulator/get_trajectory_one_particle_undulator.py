# Import libraries
import numpy as np
from numpy import sin,sqrt
from scipy.integrate import cumtrapz
from scipy.constants import c, m_e, epsilon_0, e, hbar
import h5py

# Indicate parameters for the trajectory calculation
dt_ret     = 0.01e-15                            # timestep
Tsim       = 4.0e-12                             # total simulation time
E_MeV      = 10.0                                # total electron energy (kinetic + rest) in MeV
r_beta     = 0.05e-6                             # oscillation amplitude
gamma      = E_MeV / 0.51099895                  # gamma factor
t_ret      = np.arange(0, Tsim, dt_ret)          # time array
n_e        = 5e24                                # plasma electron background density
omega_p    = sqrt(e * e * n_e/(m_e * epsilon_0)) # plasma electron frequency
omega_beta = omega_p / sqrt(2 * gamma)           # betatron frequency

# Calculate trajectory and corresponding velocities
y          = r_beta * sin(omega_beta * t_ret)    # y positions of the electron
vy         = np.zeros(np.size(y))                # prealocate y velocity of electrons

for i in range(1,len(vy)):
  vy[i]    = (y[i]-y[i-1])/dt_ret                # calculate velocity from y-positions

vy[0]    = vy[1]                                          # change the initial velocity to next time step velocity, avoid 0 value  
by       = vy / c                                         # normalized y velocity
vec_ones = np.ones(np.size(by))                           # unit vector in the size of by
bx       = np.sqrt(vec_ones-vec_ones/(gamma*gamma)-by*by) # normalized x velocity calculated from gamma and by
vx       = bx * c                                         # x velocity
x        = cumtrapz(vx, t_ret, initial=t_ret[0])          # x positions from time-integration of vx
z        = np.zeros(np.size(x))                           # z positions - zero array (no motion in z)
vz       = np.zeros(np.size(x))                           # z velocity - zero array (no motion in z)

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
with h5py.File('test_cases/one_particle_undulator/particle_trajectory.h5', 'w') as f:
    for id, arrays in data.items():
        grp = f.create_group(str(id))
        for key, values in arrays.items():
            grp.create_dataset(key, data=values)

# Print status when all done
print("HDF5 file 'particle_trajectories.h5' created successfully.")
