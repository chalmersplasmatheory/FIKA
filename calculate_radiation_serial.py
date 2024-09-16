# Import necessary libraries
from numpy import pi, sqrt
import numpy as np
from scipy.constants import c, epsilon_0, hbar, e
import h5py
import os
from scipy.interpolate import interp1d
import math
from utils import *

# Write individual particle spectra to an .h5 file
def write_spectra_to_file(output_folder, id, freq, dW_dE_dOmega, t, R0, dW_dt_dOmega, weights, weight):
  # Convert frequencies to eV energies for the output  
  freq_to_eV = 2 * pi * hbar / e 
  ene        = freq_to_eV * freq

  # Include PIC macroparticle weights in processed_data if weights are used
  if weights == True:
    processed_data = {
      "ene": ene,  
      "spectrum_ene": dW_dE_dOmega,
      "t": t-R0/c,
      "spectrum_t": dW_dt_dOmega,
      "weight":weight,
    }   
  else:
    processed_data = {
      "ene": ene,  
      "spectrum_ene": dW_dE_dOmega,
      "t": t-R0/c,
      "spectrum_t": dW_dt_dOmega,
    }  
    
  # Save the processed data under a group named after the particle ID
  with h5py.File(os.path.join(output_folder, 'individual_spectra.h5'), 'a') as outfile:
    grp = outfile.create_group(str(id))
    for data_key, data_value in processed_data.items():
      grp.create_dataset(data_key, data = data_value)

# Read, sum, and optionally PIC-weight spectra from all the particles
def read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  # Store number of particles processed to print
  curr_num_par_sum = 0 

  with h5py.File(output_folder + 'individual_spectra.h5', 'r') as file_to_sum:
    # Get the energy range of the whole spectrum from the user input 
    ene_range        = get_ene_range(E_slice_eV, file_to_sum)
    time_range       = get_time_range(t_slice_s, file_to_sum)

    # Preallocate arrays for spectrum and temporal profile
    fin_spec_ene     = np.zeros(np.size(ene_range))
    fin_spec_t       = np.zeros(np.size(time_range))

    curr_num_par_sum = 0  # Counter for processed particles

    for id in file_to_sum.keys():
      dataset          = file_to_sum[id]
      current_ene      = np.array(dataset['ene'][:])
      current_spec_ene = np.array(dataset['spectrum_ene'][:])
      current_t        = np.array(dataset['t'][:])
      current_spec_t   = np.array(dataset['spectrum_t'][:])        
        
      if weights:
        current_weight = np.array(dataset['weight'][0])
      else:
        current_weight = 1

      # Set up linear interpolators for energy and time 
      # If the value is outside the range of the interpolated array, fill_value = 0
      interpolation_ene = interp1d(current_ene, current_spec_ene * current_weight, kind='linear', bounds_error=False, fill_value=0)
      interpolation_t   = interp1d(current_t, current_spec_t * current_weight, kind='linear', bounds_error=False, fill_value=0)
      # Add the individual particle spectrum to the total final spectrum
      fin_spec_ene     += interpolation_ene(ene_range)
      fin_spec_t       += interpolation_t(time_range)
      curr_num_par_sum += 1

      if curr_num_par_sum % print_every_spectrum_sum == 0:
        print(f'Spectra of {curr_num_par_sum} particles are already summed.')
    return ene_range, fin_spec_ene, time_range, fin_spec_t

# Main function for radiation calculation
def calculate_spectrum_one_particle(charge, r, phi, theta, input_file, output_folder, weights, print_every_par_spectrum):
  # If the file spectra.h5 exists in the folder, remove it to avoid clash
  remove_existing_file(output_folder +'individual_spectra.h5')

  # Open the input file containing particle trajectories and properties
  with h5py.File(input_file, 'r') as file:
    particles_total = str(len(file))
    print(str(particles_total)+' particles will be processed.')

    # Keep track of the current particles in order to print progress
    curr_num_par = 0  

    # In case the data are not suitable for calculation, the particle will be skipped
    # Here, we count how many particles will be skipped
    particles_skipped = 0

    for id in file:
      # Extract trajectory and velocity data for the current particle
      t_ret = np.array(file[id]['time'])
      x     = np.array(file[id]['x'])
      y     = np.array(file[id]['y'])
      z     = np.array(file[id]['z'])
      vx    = np.array(file[id]['vx'])
      vy    = np.array(file[id]['vy'])
      vz    = np.array(file[id]['vz'])
      # Extract particle weight if PIC weights are being considered
      if weights == True:
        weight = np.array(file[id]['PIC_macroparticle_weight'])

      # Skip particles with insufficient data for processing
      if len(vx) < 3:
        print(f"The trajectory of particle with id {id} is too short to be processed. Skipping.")
        particles_skipped += 1
        continue

      # Calculate the time step from the original data
      dt_ret = calculate_dtime(t_ret)

      # Check if the timesteps are equally split
      # The code currently operates only with continous-like equidistant times
      continue_outer_loop = False
      for idt in range(1, len(t_ret)):
        # Check if the timestep-size differences are in the float tolerance
        if not math.isclose(t_ret[idt]-t_ret[idt-1], dt_ret):
          print('The trajectory of particle with id ' + str(id) + ' does not have equally spaced timesteps. Skipping.')
          continue_outer_loop = True # break the loop and skip the particle in case of not equal timesteps
          break
      # Check the flag to continue the outer loop
      if continue_outer_loop:
        particles_skipped += 1
        continue

      # Perform a series of calculations to determine the radiation characteristics of the particle
      beta                = calculate_beta(vx, vy, vz)
      gamma               = 1/sqrt(1-beta[0]*beta[0]-beta[1]*beta[1]-beta[2]*beta[2])
      beta_dot            = calculate_beta_dot(beta, dt_ret )
      pos_vec             = calculate_pos_vec(r, phi, theta)  
      R0                  = calculate_R0(pos_vec)
      n                   = calculate_n(pos_vec, R0) 
      R                   = calculate_R(x, y, z, r, phi, theta)  
      t                   = calculate_t(t_ret, R)  
      E_x, E_y, E_z       = calculate_E(charge, R, n, beta, beta_dot)
      signal_x            = R * E_x
      signal_y            = R * E_y
      signal_z            = R * E_z

      # Interpolate the time signal to a better resolution
      # The timestep is given by the particle energy and original timestep
      # This captures the highest photon frequencies
      # See for instance M. Pardal et al.,Computer Physics Communications, 285, 2023, 108634 for explanation
      dt_inter            = 0.5 * dt_ret / max(gamma)**2
      
      Nt                  = calculate_number_frequency_steps(t, dt_inter)
      t_inter             = np.linspace(min(t), max(t), Nt)
      signal_inter_x      = np.interp(t_inter, t, signal_x)
      signal_inter_y      = np.interp(t_inter, t, signal_y)
      signal_inter_z      = np.interp(t_inter, t, signal_z)
      signal_squared      = np.square(signal_x) + np.square(signal_y) + np.square(signal_z)      
      dW_dt_dOmega        = c * epsilon_0 * signal_squared 
      freq, fft_squared_x = calculateFFTsquared(dt_inter, signal_inter_x)
      freq, fft_squared_y = calculateFFTsquared(dt_inter, signal_inter_y)
      freq, fft_squared_z = calculateFFTsquared(dt_inter, signal_inter_z)
      fft_squared         = fft_squared_x + fft_squared_y + fft_squared_z
      dW_dE_dOmega        = calculate_spectrum(fft_squared) 

      # Save the calculated spectra for the particle to an HDF5 file
      if weights == True:
        write_spectra_to_file(output_folder, id, freq, dW_dE_dOmega, t, R0, dW_dt_dOmega, weights, weight)
      else:
        write_spectra_to_file(output_folder, id, freq, dW_dE_dOmega, t, R0, dW_dt_dOmega, weights, False)

      # Print progress at regular intervals
      print_calculation_progress(curr_num_par, print_every_par_spectrum)
      curr_num_par += 1 
    return particles_total, particles_skipped

# Function to aggregate spectra from all particles and save the final results
# Summarizes the collective radiation impact from all particles
def calculate_spectrum_all_particles(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  remove_existing_file(output_folder + 'final_spectrum.h5')
  ene_range, fin_spec_ene, time_range, fin_spec_t = read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum)
  write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, ene_range, fin_spec_ene)