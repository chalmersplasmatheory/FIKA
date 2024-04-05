# Import necessary libraries
from numpy import pi, sin, cos, sqrt
import numpy as np
from scipy.constants import c, epsilon_0, hbar, e
import h5py
import os

# Calculate the time step based on an array of time values
def calculate_dtime(time):
  dtime      = time[1] - time[0] 
  return dtime

# Calculate the dimensionless velocity (beta) components of a particle
def calculate_beta(vx,vy,vz):
  beta_x   = np.array(vx / c)
  beta_y   = np.array(vy / c)
  beta_z   = np.array(vz / c)
  beta     = np.array([beta_x, beta_y, beta_z])  # Combine into a single array
  return beta

# Calculate the derivative of beta with respect to time (acceleration)
def calculate_beta_dot(beta, dt_ret):
  beta_x     = beta[0, :]
  beta_y     = beta[1, :]
  beta_z     = beta[2, :]
  beta_x_dot = (beta_x[1:-1] - beta_x[0:-2]) / dt_ret  
  beta_x_dot = np.append(beta_x_dot, [beta_x_dot[-1], beta_x_dot[-1]])           
  beta_y_dot = (beta_y[1:-1] - beta_y[0:-2]) / dt_ret       
  beta_y_dot = np.append(beta_y_dot, [beta_y_dot[-1], beta_y_dot[-1]])               
  beta_z_dot = (beta_z[1:-1] - beta_z[0:-2]) / dt_ret       
  beta_z_dot = np.append(beta_z_dot, [beta_z_dot[-1], beta_z_dot[-1]]) 
  beta_dot   = np.array([beta_x_dot, beta_y_dot, beta_z_dot]) 
  return beta_dot

# Calculate position vector (of the observer from the axis origin) from spherical coordinates 
def calculate_pos_vec(r, phi, theta):  
  pos_vec  = np.array([r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)])
  return pos_vec

# Calculate the magnitude of the position vector (of the observer from the axis origin)
def calculate_R0(pos_vec):  
  R0       = np.linalg.norm(pos_vec)    
  return R0

# Calculate the unit vector in the direction of observation
def calculate_n(pos_vec, R0):                                  
  n        = np.array(pos_vec / R0)
  return n

# Calculate distance of the particle from the observer
def calculate_R(x, y, z, r, phi, theta):                             
  pos_vec = calculate_pos_vec(r, phi, theta)
  R       = sqrt((pos_vec[0] - x)**2 + (pos_vec[1] - y)**2 + (pos_vec[2] - z)**2)
  return R

# Calculate observer's time t
def calculate_t(t_ret,R):                                  
  t        =  t_ret + R / c 
  return t

# Calculate electric field at the observer's position
def calculate_E(charge,R,n,beta,beta_dot):
  beta_x      = beta[0, :]
  beta_y      = beta[1, :]
  beta_z      = beta[2, :]                               
  Econst      = charge / (4 * pi * epsilon_0 * c) 
  n_vecs      = np.tile(n, (np.size(beta_x), 1)).T
  n_min_b     = n_vecs - beta
  n_beta      = n[0] * beta_x + n[1] * beta_y + n[2] * beta_z; 
  numerator   = np.cross(n.T, np.cross(n_min_b.T, beta_dot.T)).T
  denominator = (np.power(1 - n_beta, 3) * R)
  E_x         = Econst * numerator[0, :] / denominator
  E_y         = Econst * numerator[1, :] / denominator
  E_z         = Econst * numerator[2, :] / denominator
  return E_x, E_y, E_z 

# Determine the number of frequency steps from the maximum photon energy given by a user
# This function calculates the sampling frequency necessary to resolve the spectrum up to a maximum energy defined by the user
def calculate_number_frequency_steps(t, E_radMax_eV):
  E_radMax      = E_radMax_eV * e        # Convert energy from eV to Joules
  omega_rad_max = E_radMax / hbar        # Convert energy to angular frequency
  freq_rad_max  = omega_rad_max / 2 / pi # Convert angular frequency to frequency
  t_inter       = 1/ (20 * freq_rad_max) # Define time interval for 20 samples per period of the highest frequency
  Nt            = int(round((max(t)-min(t))/t_inter))
  return Nt

# Perform Fast Fourier Transform on a signal and return squared magnitude of its frequency components and corresponding frequencies
def calculateFFTsquared(dt_inter, signal_inter):
  N                           = len(signal_inter) # Number of points in the signal
  fft_output                  = np.fft.fft(signal_inter)
  frequencies_whole_signal    = np.fft.fftfreq(N, dt_inter)
  # Normalize FFT output to make its amplitude independent of the number of points
  normalized_fft_whole_signal = fft_output / N
  # Return the first half of the frequency spectrum and squared magnitude (for power spectrum)
  # Output only 1/2 the frequency spectrum, it is sufficient for real-valued signals because the FFT output is symmetric 
  return frequencies_whole_signal[:N // 2], np.square(np.abs(normalized_fft_whole_signal)[:N // 2])

# Calculate spectral intensity dI/dE/dOmega in units J/eV/sr
def calculate_spectrum(fft_squared):
  # c * epsilon_0 / pi comes from the radiation formula
  # hbar / e is for transforming frequency to eV
  Iconst = (hbar / e) * (c * epsilon_0 / pi)
  # Apply Parseval's theorem to ensure energy conservation in FFT, ...
  # ... doubling the intensity to account for only using half the spectrum
  return Iconst * fft_squared * 2

# Write individual particle spectra to an .h5 file
def write_spectra_to_file(output_folder, id, freq, dI_dE_dOmega, t, R0, dI_dt_dOmega, weights, weight):
  # Include PIC macroparticle weights in processed_data if weights are used
  if weights ==True:
    processed_data = {
      "freq": freq,  
      "spectrum_freq": dI_dE_dOmega,
      "t": t-R0/c,
      "spectrum_t": dI_dt_dOmega,
      "weight":weight,
    }   
  else:
    processed_data = {
      "freq": freq,  
      "spectrum_freq": dI_dE_dOmega,
      "t": t-R0/c,
      "spectrum_t": dI_dt_dOmega,
    }  
    
  # Save the processed data under a group named after the particle ID
  with h5py.File(os.path.join(output_folder, 'spectra.h5'), 'a') as outfile:
    grp = outfile.create_group(str(id))
    for data_key, data_value in processed_data.items():
      grp.create_dataset(data_key, data = data_value)

# Print progress of calculations for every "print_freq" particles processed
def print_calculation_progress(curr_num_par,print_freq):
  if curr_num_par % print_freq == 0  and curr_num_par != 0:
    print(str(curr_num_par)+' particles are already processed.')

# Determine the time range for the summation of all the spectra
# This function finds the overall time interval across all particles to aggregate the spectra consistently
def get_time_range(t_slice_s, file_to_sum): 
  t_id_index = 0
  for par_id in file_to_sum.keys():
    # Read the observer times for the current particle
    current_t     = np.array(file_to_sum[par_id]['t'][:])
    # Maximum & minimum time of the observer signal for the current particle
    current_max_t = np.max(current_t)
    current_min_t = np.min(current_t)
    # Set the minimum and maximum time to the values of the 1st particle, then update if needed
    if t_id_index == 0 or current_max_t > max_t:
      max_t = current_max_t             
    if t_id_index == 0 or current_min_t < min_t:
      min_t = current_min_t  
    t_id_index += 1    
  
  # The time range of the whole signal is from min to max observer time from all the particles 
  num_t      = int(round((max_t - min_t) / t_slice_s))
  # Create an array of time steps within this range
  time_range = np.linspace(min_t, max_t, num = num_t)
  return(time_range)

# Determine the frequency range for the summation of all the spectra
# Similar to get_time_range, this function finds the overall frequency range for aggregating spectral data
def get_freq_range(E_slice_eV, file_to_sum):
  freq_id_index = 0
  for par_id in file_to_sum.keys():
    # Read the frequencies for the current particle
    current_freq     = np.array(file_to_sum[par_id]['freq'][:])
    # Maximum frequency for the current particle
    current_max_freq = np.max(current_freq)
    # Set the  maximum frequency to the values of the 1st particle, then update if needed
    if freq_id_index == 0 or current_max_freq > max_freq:
      max_freq = current_max_freq
    freq_id_index += 1   

  # The user energy input is converted to the unit of 1/s
  freq_un        = e / (2*pi*hbar) * E_slice_eV
  # The upper frequency is determined as max from all the particles 
  num_freq       = int(round(max_freq/freq_un))
  freq_range     = np.linspace(0, max_freq, num = num_freq)
  return freq_range  

# Read, sum, and optionally PIC-weight spectra from all the particles
def read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  # Store number of particles processed to print
  curr_num_par_sum = 0 

  with h5py.File(output_folder + 'spectra.h5', 'r') as file_to_sum:
    # Get the frequency range of the whole spectrum from the user input 
    freq_range    = get_freq_range(E_slice_eV, file_to_sum)
    time_range    = get_time_range(t_slice_s, file_to_sum)

    # Preallocate arrays for spectrum and temporal profile
    fin_spec_freq  = np.zeros(np.size(freq_range))
    fin_spec_t     = np.zeros(np.size(time_range))

    for id in file_to_sum.keys():
      dataset        = file_to_sum[id]
      current_freq   = np.array(dataset['freq'][:])
      current_spec_freq   = np.array(dataset['spectrum_freq'][:])
      current_t      = np.array(dataset['t'][:])
      current_spec_t = np.array(dataset['spectrum_t'][:])
      if weights:
        current_weight = np.array(dataset['weight'][0])

      # Sum the spectral data, accounting for the weight of each particle if applicable
      # Each particle has a different frequency sampling
      # The intensity value is assigned to the closest frequency value of the final spectrum sampling
      for f_num in range(0, len(current_freq)):
        index_of_closest_freq = (np.abs(freq_range - current_freq[f_num])).argmin()
        # If weights are included, spectrum from each macroparticle is multiplied by the real number of particles
        if weights == True:
          current_spec_to_sum = current_spec_freq[f_num] * current_weight
        else:
          current_spec_to_sum = current_spec_freq[f_num]            
        fin_spec_freq[index_of_closest_freq] += current_spec_to_sum

      # Sum the spectral data for observer time signal, similarly to the frequency above
      for t_num in range(0, len(current_t)):
        index_of_closest_t              = (np.abs(time_range - current_t[t_num])).argmin()
        if weights == True:
          current_spec_t_to_sum = current_spec_t[t_num] * current_weight
        else:
          current_spec_t_to_sum = current_spec_t[t_num]            
        fin_spec_t[index_of_closest_t] += current_spec_t_to_sum
      
      # Print progress at intervals  
      if curr_num_par_sum % print_every_spectrum_sum == 0 and curr_num_par_sum != 0:
        print('Spectrum from ' + str(curr_num_par_sum)+' particles is already summed.')
      curr_num_par_sum += 1
  return freq_range, fin_spec_freq, time_range, fin_spec_t

# For writing the summation of all particle spectra into the final output file
# Saves the aggregated spectra and temporal profiles to a final HDF5 file
def write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, freq_range, fin_spec_freq):
  with h5py.File(output_folder + 'final_spectrum.h5', 'w') as final_file_to_write:
    final_file_to_write.create_dataset('freq', data = freq_range )
    final_file_to_write.create_dataset('spectrum_freq', data = fin_spec_freq)   
    final_file_to_write.create_dataset('t', data = time_range )
    final_file_to_write.create_dataset('spectrum_t', data = fin_spec_t)     

# For removing old output files before starting a new calculation
# Checks if a file exists at a specified path and removes it to prevent data overlap
def remove_existing_file(file_path):    
  file_exists = os.path.exists(file_path)
  if file_exists:
    os.remove(file_path)
    print('File ' + str(file_path) + ' already existed. The old file was removed.')

# Main function for radiation calculation
def calculate_spectrum_one_particle(charge,E_radMax_eV, r, phi, theta, input_file, output_folder, weights, print_every_par_spectrum):

  # If the file spectra.h5 exists in the folder, remove it to avoid clash
  remove_existing_file(output_folder +'spectra.h5')

  # Open the input file containing particle trajectories and properties
  with h5py.File(input_file, 'r') as file:
    print(str(len(file))+' particles will be processed.')

    # Keep track of the current particles in order to print progress
    curr_num_par = 0  

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
        continue

      # Perform a series of calculations to determine the radiation characteristics of the particle
      dt_ret              = calculate_dtime(t_ret)
      beta                = calculate_beta(vx, vy, vz)
      beta_dot            = calculate_beta_dot(beta, dt_ret )
      pos_vec             = calculate_pos_vec(r, phi, theta)  
      R0                  = calculate_R0(pos_vec)
      n                   = calculate_n(pos_vec, R0) 
      R                   = calculate_R(x,y,z,r,phi, theta)  
      t                   = calculate_t(t_ret,R)  
      E_x, E_y, E_z       = calculate_E(charge, R, n, beta, beta_dot)
      signal_x            = R * E_x
      signal_y            = R * E_y
      signal_z            = R * E_z
      Nt                  = calculate_number_frequency_steps(t, E_radMax_eV)
      t_inter             = np.linspace(min(t),max(t),Nt)
      signal_inter_x      = np.interp(t_inter, t, signal_x)
      signal_inter_y      = np.interp(t_inter, t, signal_y)
      signal_inter_z      = np.interp(t_inter, t, signal_z)
      signal_squared      = np.square(signal_x) + np.square(signal_y) + np.square(signal_z)      
      dI_dt_dOmega        = c * epsilon_0 * signal_squared 
      dt_inter            = calculate_dtime(t_inter)
      freq, fft_squared_x = calculateFFTsquared(dt_inter, signal_inter_x)
      freq, fft_squared_y = calculateFFTsquared(dt_inter, signal_inter_y)
      freq, fft_squared_z = calculateFFTsquared(dt_inter, signal_inter_z)
      fft_squared         = fft_squared_x + fft_squared_y + fft_squared_z
      dI_dE_dOmega        = calculate_spectrum(fft_squared) 

      # Save the calculated spectra for the particle to an HDF5 file
      if weights == True:
        write_spectra_to_file(output_folder, id, freq, dI_dE_dOmega, t, R0, dI_dt_dOmega, weights, weight)
      else:
        write_spectra_to_file(output_folder, id, freq, dI_dE_dOmega, t, R0, dI_dt_dOmega, weights, False)

      # Print progress at regular intervals
      print_calculation_progress(curr_num_par, print_every_par_spectrum)
      curr_num_par += 1 

# Function to aggregate spectra from all particles and save the final results
# Summarizes the collective radiation impact from all particles
def calculate_spectrum_all_particles(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  freq_range, fin_spec_freq, time_range, fin_spec_t = read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum)
  remove_existing_file(output_folder + 'final_spectrum.h5')
  write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, freq_range, fin_spec_freq)
