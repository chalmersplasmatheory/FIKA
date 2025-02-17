# Import necessary library
from numpy import pi, sin, cos, sqrt
import numpy as np
from scipy.constants import c, epsilon_0, hbar, e
import h5py
import os

# Print progress of calculations for every "print_freq" particles processed
def print_calculation_progress(curr_num_par,print_freq):
  if curr_num_par % print_freq == 0  and curr_num_par != 0:
    print(str(curr_num_par)+' particles are already processed.')

# For removing old output files before starting a new calculation
# Checks if a file exists at a specified path and removes it to prevent data overlap
def remove_existing_file(file_path):    
  file_exists = os.path.exists(file_path)
  if file_exists:
    os.remove(file_path)
    print('File ' + str(file_path) + ' already existed. The old file was removed and replaced by a new one.')

# Calculate the time step based on an array of time values
def calculate_dtime(time):
  dtime      = time[1] - time[0] 
  return dtime

# Calculate the dimensionless velocity (beta) components of a particle
def calculate_beta(vx, vy, vz):
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
  beta_x_dot = np.append(beta_x_dot, [beta_x_dot[-1], beta_x_dot[-1]]) # append 2x with the last element to fit the size        
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
def calculate_E(charge, R, n, beta, beta_dot):
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

# Determine the number of frequency steps from the interpolation timestep
def calculate_number_frequency_steps(t, dt_inter):
  Nt            = int(round((max(t)-min(t))/dt_inter))
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

# Calculate spectral intensity dW/dE/dOmega in units J/eV/sr
def calculate_spectrum(fft_squared):
  # c * epsilon_0 / pi comes from the radiation formula
  Wconst = (c * epsilon_0 / pi)
  # Apply Parseval's theorem to ensure energy conservation in FFT, ...
  # ... doubling the intensity to account for only using half the spectrum
  return Wconst * fft_squared * 2

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

# Determine the energy range for the summation of all the spectra
# Similar to get_time_range, this function finds the overall energy range for aggregating spectral data
def get_ene_range(E_slice_eV, file_to_sum):
  ene_id_index = 0
  for par_id in file_to_sum.keys():
    # Read the energies for the current particle
    current_ene     = np.array(file_to_sum[par_id]['ene'][:])
    # Maximum energy for the current particle
    current_max_ene = np.max(current_ene)
    # Set the maximum energy to the values of the 1st particle, then update if needed
    if ene_id_index == 0 or current_max_ene > max_ene:
      max_ene = current_max_ene
    ene_id_index += 1   

  # The upper requency is determined as max from all the particles 
  num_ene       = int(round(max_ene/E_slice_eV))
  ene_range     = np.linspace(0, max_ene, num = num_ene)
  return ene_range  

# For writing the summation of all particle spectra into the final output file
# Saves the aggregated spectra and temporal profiles to a final HDF5 file
def write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, ene_range, fin_spec_ene):
  with h5py.File(output_folder + 'final_spectrum.h5', 'w') as final_file_to_write:
    final_file_to_write.create_dataset('ene', data = ene_range)
    final_file_to_write.create_dataset('spectrum_ene', data = fin_spec_ene)   
    final_file_to_write.create_dataset('t', data = time_range)
    final_file_to_write.create_dataset('spectrum_t', data = fin_spec_t)     
