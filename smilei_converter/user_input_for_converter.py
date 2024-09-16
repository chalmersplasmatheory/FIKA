# THIS IS AN EXAMPLE INPUT FILE FOR SMILEI CONVERTER

# Set the reference angular frequency from your Smilei simulation.
# Set this value in the SI units rad/s (not Smilei units !!!).
# Example: omega_r_smilei_SI = 2354564459136066.5 is for frequency 2354564459136066.5 rad/s, corresponding to 0.8 micron wavelength.
omega_r_smilei_SI      = 2354564459136066.5

# Set the timestep used in your Smilei simulation (i.e. the timestep of the PIC loop).
# It is necessary that you convert this value to the SI units !!!
# Example: dt_smilei_SI = 8.539240e-17 is for timestep 8.539240e-17 s.
dt_smilei_SI           =  8.539240e-17

# Set the full path to the Smilei .h5 file containing particle tracking data, including the file name.
# It has to be the file that is already sorted by postprocessing with happi. Unsorted versions cannot be processed !!!
smilei_file_to_convert = 'INTENSE/test_cases/PIC_trajectories/TrackParticles_electron.h5'

# Set the full path of the output file .h5 where the data for radiation calculation will be stored.
converted_file_name    = 'INTENSE/test_cases/PIC_trajectories/test_particle_set.h5'

# Progress update frequency. After new set of "print_every" particles is processed, the update will be printed into a standard output.
print_every            = 100
