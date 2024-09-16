# Import necessary libraries
from numpy import pi, sqrt
import numpy as np
from scipy.constants import c, epsilon_0, hbar, e
import h5py
import os
from scipy.interpolate import interp1d
import math
from utils import *
from mpi4py import MPI
import sys

# NOTE: Ensure that H5PY is compiled with the parallel version for MPI compatibility

# Initialize MPI
# Setting up a communicator (a group of processes that can communicate with each other)
comm              = MPI.COMM_WORLD
# Determine the rank of the current process 
rank              = comm.Get_rank()
# Find out the total number of processes that are part of the communicator comm
size              = comm.Get_size()
# Rank of the root process
root_process_rank = 0

# Read, sum, and optionally PIC-weight spectra from all the particles
def read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  # Open the HDF5 file containing individual spectra with MPI support
  with h5py.File(output_folder + 'individual_spectra.h5', 'r', driver='mpio', comm=comm) as file_to_sum:
    # Get energy and time ranges from the file
    ene_range          = get_ene_range(E_slice_eV, file_to_sum)
    time_range         = get_time_range(t_slice_s, file_to_sum)

    # Initialize local arrays for summing spectra
    local_spec_ene     = np.zeros(np.size(ene_range))
    local_spec_t       = np.zeros(np.size(time_range))
    local_par_num      = 0

    # Distribute keys among processes
    keys               = list(file_to_sum.keys())
    local_keys         = keys[rank::size]  # Divide keys among processes
    local_par_num_curr = 0
    
    # Only the root process (rank 0) initializes the print point variable
    if rank == 0:
      print_point = print_every_spectrum_sum

    # Loop through assigned keys to read and sum spectra
    for key in local_keys:
      dataset          = file_to_sum[key]
      current_ene      = np.array(dataset['ene'][:])
      current_spec_ene = np.array(dataset['spectrum_ene'][:])
      current_t        = np.array(dataset['t'][:])
      current_spec_t   = np.array(dataset['spectrum_t'][:])
      local_par_num   += 1

      # Apply weights if specified
      if weights:
        current_weight = np.array(dataset['weight'][0])
      else:
        current_weight = 1

      # Interpolate spectra to the defined energy and time ranges
      # Each particle's spectra are recorded over different energy and time ranges.
      # To sum these spectra coherently, we need to interpolate them onto a common set of ranges.

      # Interpolating energy spectrum:
      # Create an interpolation function for the energy spectrum of the current particle.
      # - `current_ene`: The energy points of the current particle's spectrum.
      # - `current_spec_ene`: The spectrum values at these energy points.
      # - `current_weight`: The weight applied to the current spectrum (either 1 or a specified weight).
      # - `kind='linear'`: Use linear interpolation.
      # - `bounds_error=False`: Prevent errors for values outside the interpolation range.
      # - `fill_value=0`: Fill with zeros for out-of-range values.
      interpolation_ene   = interp1d(current_ene, current_spec_ene * current_weight, kind='linear', bounds_error=False, fill_value=0)
      interpolation_t     = interp1d(current_t, current_spec_t * current_weight, kind='linear', bounds_error=False, fill_value=0)

      # Apply the interpolation functions to the common energy and time ranges:
      # The interpolated values are added to the local MPI-process spectrum sums.
      local_spec_ene     += interpolation_ene(ene_range)
      local_spec_t       += interpolation_t(time_range)

      # Increment the count of processed particles for the current process:
      local_par_num_curr += 1
      local_par_num       = np.array([local_par_num_curr], dtype=np.float64)
      # Prepare an array to hold the global particle count after reduction:
      global_par_num      = np.zeros(1)

      # Collect the local particle counts from all processes and sums them up.
      comm.Reduce(np.array([local_par_num]), global_par_num, op=MPI.SUM, root=0)
 
      if rank == 0:
        if global_par_num[0] > print_point:
          print(f'Spectra of {int(print_point)} particles are already summed.',flush = True)
          print_point += print_every_spectrum_sum
  # Ensure that all processes are finished before the final summation
  comm.Barrier()

  # Prealocate the arrays for final spectrum summation
  global_spec_ene = np.zeros_like(local_spec_ene)
  global_spec_t    = np.zeros_like(local_spec_t)
  global_par_num   = np.array([0.], dtype=np.float64)

  #Sum spectra from all the processes
  comm.Reduce(local_spec_ene, global_spec_ene, op=MPI.SUM, root=0)
  comm.Reduce(local_spec_t, global_spec_t, op=MPI.SUM, root=0)
  comm.Reduce(local_par_num, global_par_num, op=MPI.SUM, root=0)
  comm.Reduce(np.array([local_par_num]), global_par_num, op=MPI.SUM, root=0)

  # On the root process, print the total number of particles summed
  if rank == 0:
    print("Calculation finished. Spectra from", int(global_par_num[0]), "particles were summed.")
    return ene_range, global_spec_ene, time_range, global_spec_t
  else:
    # Explicitly return nothing on non-root processes
    return None, None, None, None

def write_data_sequentially(key, processed_data, output_folder):
  # This function will be called by process 0 to write data sequentially
  with h5py.File(os.path.join(output_folder, 'individual_spectra.h5'), 'a') as outfile:
    # Check if the key already exists in the file
    if key not in outfile:
      # Create a new group for the key
      grp = outfile.create_group(key)
      # Write spectrum for each key
      for data_key, data_value in processed_data.items():
        grp.create_dataset(data_key, data = data_value)

# Count the number of particles in the input file and create a new HDF5 file for writing.
def count_particles_and_create_file_to_write(input_file,output_folder):  
  # Only the root process executes this block
  if rank == 0:
    with h5py.File(input_file, 'r') as file:
      par_num = len(file)
      file_path = os.path.join(output_folder, 'individual_spectra.h5')
    try:
        # Attempt to remove the existing HDF5 file if it exists
        os.remove(file_path)
        print(f"The original {file_path} has been deleted, creating a new one.")
    except FileNotFoundError:
        print(f"The file {file_path} was not found, creating a new one.")
    return par_num

# Main function for radiation calculation
def calculate_spectrum_one_particle(charge, r, phi, theta, input_file, output_folder, weights, print_every_par_spectrum):
  local_skipped = 0
  par_num = count_particles_and_create_file_to_write(input_file,output_folder)
  if rank == 0:
    print(str(par_num) + ' particles will be processed.', flush=True)
  comm.Barrier()
  keys = None
  if rank == 0:
    with h5py.File(os.path.join(input_file), 'r') as infile:
      keys = sorted(list(infile.keys()))
  comm.Barrier()
  keys                 = comm.bcast(keys, root=0)  # Broadcasting the keys to all processes
  par_num_to_count     = 0
  tot_par_num_to_count = 0
  comm.Barrier()
  for i in range(0, len(keys), size):
    chunk_keys = keys[i:i + size]   
    # Each process reads and processes its assigned key in the current chunk
    if rank < len(chunk_keys):
      key_to_process = chunk_keys[rank]
    else:       
      key_to_process = None

    with h5py.File(os.path.join(input_file), 'r', driver='mpio', comm=comm) as file:
      if key_to_process == None:
        processed_data=None
      else:
        key            = key_to_process
        particle_id    = key
        particle_group = file[particle_id]
        t_ret          = np.array(particle_group['time'])
        x, y, z        = np.array(particle_group['x']), np.array(particle_group['y']), np.array(particle_group['z'])
        vx, vy, vz     = np.array(particle_group['vx']), np.array(particle_group['vy']), np.array(particle_group['vz'])
        if weights == True:
          weight         = np.array(particle_group['PIC_macroparticle_weight'])

        skip_particle = False
        # Check for sufficient trajectory data
        if len(vx) < 3:
          print(f"The trajectory of particle with id {particle_id} is too short to be processed. Skipping on rank {rank}.")
          skip_particle = True  # Local flag to skip processing         

        if not skip_particle:
          # Calculate the time step from the original data
          dt_ret = calculate_dtime(t_ret)
          # Check if the timesteps are equally split
          for idt in range(1, len(t_ret)):
            if not math.isclose(t_ret[idt] - t_ret[idt - 1], dt_ret):
              print(f'The trajectory of particle with id {particle_id} does not have equally spaced timesteps. Skipping on rank {rank}.')
              skip_particle = True  # Local flag to skip processing

        if skip_particle:
          processed_data = None
          local_skipped += 1 
        else:
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
              "weight": weight,
            }   
          else:
            processed_data = {
              "ene": ene,  
              "spectrum_ene": dW_dE_dOmega,
              "t": t-R0/c,
              "spectrum_t": dW_dt_dOmega,
            }  
        
      all_processed_data = comm.gather(processed_data, root=0) 
      if rank == 0:        
          par_num_to_count += len(all_processed_data)

          for j, key in enumerate(chunk_keys):
            if not all_processed_data[j] == None:
              write_data_sequentially(key, all_processed_data[j], output_folder)
          if par_num_to_count >= print_every_par_spectrum: 
            tot_par_num_to_count += par_num_to_count
            print("Spectrum from",tot_par_num_to_count,"particles has been already calculated.")
            par_num_to_count= 0   
    comm.Barrier()  
  comm.Barrier()    
  particles_skipped = comm.reduce(local_skipped, op = MPI.SUM, root = 0)
  if rank == 0:
    print('Calculation of the spectra of individual particles is finished. '+ str(particles_skipped)+ ' particles of the original number '+ str(len(keys))+ ' were skipped in the calculation.', flush=True)
    print('Radiation spectra from '+ str(len(keys) - particles_skipped) +' particles was calculated.')

# Function to aggregate spectra from all particles and save the final results
# Summarizes the collective radiation impact from all particles
def calculate_spectrum_all_particles(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  ene_range, fin_spec_ene, time_range, fin_spec_t = read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum)
  if rank == 0:
    remove_existing_file(output_folder + 'final_spectrum.h5')
    write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, ene_range, fin_spec_ene)

# Write the final spectrum data to an HDF5 file
def write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, ene_range, fin_spec_ene):
  with h5py.File(os.path.join(output_folder, 'final_spectrum.h5'), 'w') as final_file_to_write:
    final_file_to_write.create_dataset('ene', data = ene_range)
    final_file_to_write.create_dataset('spectrum_ene', data = fin_spec_ene)
    final_file_to_write.create_dataset('t', data = time_range)
    final_file_to_write.create_dataset('spectrum_t', data=fin_spec_t)

# Read spectra from all the particles from file and sum them
def calculate_spectrum_all_particles(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum):
  if rank == 0:
    print('Summation of the spectra starts.')
  ene_range, fin_spec_ene, time_range, fin_spec_t = read_and_sum_spectra(output_folder, E_slice_eV, t_slice_s, weights, print_every_spectrum_sum)
  if rank == 0:
    write_final_spectrum_to_file(output_folder, time_range, fin_spec_t, ene_range, fin_spec_ene)
