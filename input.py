# THIS IS AN EXAMPLE OF THE INPUT FILE
# FOR DETAILED INFORMATION ON THE INPUT VARIABLES FOLLOW README.md
# Import constants from Python libraries
from numpy import pi 
from scipy.constants import e

input_file                = 'test_cases/PIC_trajectories/test_particle_set.h5' # Input file
output_folder             = 'test_cases/PIC_trajectories/no_weights_parallel/' # Output folder, include / at the end  
PIC_macroparticle_weights = False                                              # Include PIC macroparticle weights or not
parallel                  = True                                               # Serial or parallel computing
charge                    = -e                                                 # Particle charge                                                                # Max photon energy for Fourier sampling, the sampling is 20x this value
r                         = 1                                                  # Radial coordinate of the observer in m (spherical coordinates)
phi                       = 0                                                  # Azimuthal angle of observer in rad (spherical coordinates)
theta                     = pi/2                                               # Polar angle of observer in rad (spherical coordinates)
print_every_par_spectrum  = 1000                                               # Progress update frequency for individual particle spectra calculated
sum_spectra               = True                                               # After spectra fo individual particles are calculated, the summation starts or not

# INPUT FOR THE SPECTRUM SUMMATION
# NOT NEEDED IF sum_spec = False
E_slice_eV                = 10                                                 # Energy bin size in eV for final energy histogram
t_slice_s                 = 1e-19                                              # Time bin size in s for final temporal profile 
print_every_spectrum_sum  = 1000                                               # Progress update frequency for spectrum summation
