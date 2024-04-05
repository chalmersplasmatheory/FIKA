# THIS IS AN EXAMPLE OF THE INPUT FILE
# FOR DETAILED INFORMATION ON THE INPUT VARIABLES FOLLOW README.md

# Import constants from Python libraries
from numpy import pi 
from scipy.constants import e

# INPUT VARIABLES FOR SPECTRA CALCULATIONS
input_file                = 'test_case_1/test_particle_set_1.h5'         # Input file
output_folder             = 'test_case_1/with_macroparticle_weights/'    # Where to convert particle-in-cell macroparticle spectrum to the real particle spectrum.  
PIC_macroparticle_weights = True                                         # Include PIC macroparticle weights or not
parallel                  = 'Disabled'                                   # Serial or parallel computing
charge                    = -e                                           # Particle charge 
E_radMax_eV               = 1000                                         # Max photon energy for Fourier sampling, the sampling is 20x this value
r                         = 10000                                        # Radial coordinate of the observer in m (spherical coordinates)
phi                       = 0                                            # Azimuthal angle of observer in rad (spherical coordinates)
theta                     = pi/2                                         # Polar angle of observer in rad (spherical coordinates)
print_every_par_spectrum  = 200                                          # Progress update frequency for individual particle spectra calculated
sum_spectra               = True                                         # After spectra fo individual particles are calculated, the summation starts or not

# INPUT FOR THE SPECTRUM SUMMATION
# NOT NEEDED IF sum_spec = False
E_slice_eV                = 2                                            # Energy bin size in eV for final energy histogram
t_slice_s                 = 1e-16                                        # Time bin size in s for final temporal profile 
print_every_spectrum_sum  = 200                                          # Progress update frequency for spectrum summation
