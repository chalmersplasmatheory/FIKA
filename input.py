# THIS IS AN EXAMPLE OF THE INPUT FILE
# FOR DETAILED INFORMATION ON THE INPUT VARIABLES FOLLOW README.md
# Import constants from Python libraries
from numpy import pi 
from scipy.constants import e
import os

# Create a string representing the closest common directory
base_folder = os.path.dirname(os.path.abspath(__file__))

input_file                = base_folder + '/test_cases/one_particle_circular_orbit/particle_trajectory.h5' # Input file
output_folder             = base_folder + '/test_cases/one_particle_circular_orbit/' # Output folder, include / at the end  
PIC_macroparticle_weights = False                                              # Include PIC macroparticle weights or not
parallel                  = False                                              # Serial or parallel computing
charge                    = -e                                                 # Particle charge                                                                # Max photon energy for Fourier sampling, the sampling is 20x this value
r                         = 1                                                  # Radial coordinate of the observer in m (spherical coordinates)
phi                       = pi/4                                                  # Azimuthal angle of observer in rad (spherical coordinates)
theta                     = pi/2                                               # Polar angle of observer in rad (spherical coordinates)
print_every_par_spectrum  = 1000                                               # Progress update frequency for individual particle spectra calculated
sum_spectra               = True                                               # After spectra fo individual particles are calculated, the summation starts or not

# INPUT FOR THE SPECTRUM SUMMATION
# NOT NEEDED IF sum_spec = False
E_slice_eV                = 0.1                                                # Energy bin size in eV for final energy histogram
t_slice_s                 = 1e-19                                              # Time bin size in s for final temporal profile 
print_every_spectrum_sum  = 1000                                               # Progress update frequency for spectrum summation
zero_padding_res_factor   = 1                                                # Set to a value between in the interval (0,1] to control the resolution
