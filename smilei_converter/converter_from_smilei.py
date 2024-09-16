# Import necessary libraries
import h5py
import numpy as np
from numpy import pi
from scipy.constants import c, epsilon_0, m_e, e

def convert_smilei_file(omega_r_smilei_SI, dt_smilei_SI, smilei_file_to_convert, converted_file_name, print_every):
  #   Convert particle data from a Smilei simulation into SI units and store it in a new HDF5 file.
  #   Parameters:
  #      omega_r_smilei_SI (float): Angular frequency reference in Smilei, in SI units (rad/s).
  #      dt_smilei_SI (float): Time step of the simulation in SI units (s).
  #      smilei_file_to_convert (str): Path to the input Smilei HDF5 file.
  #      converted_file_name (str): Path for the output HDF5 file with converted data.
  #      print_every (int): Interval at which to print progress updates (based on particle number).

  def get_conversion_factors(omega_r_smilei_SI, dt_smilei_SI):
    # Calculate SI unit factors based on the simulation parameters
    time_SI_fac    = dt_smilei_SI
    lambda_r_SI    = 2.0 * pi * c / omega_r_smilei_SI                  # convert the reference wavelength to SI
    lamba_r_smilei = 2.0 * pi                                          # reference wavelength in Smilei units
    Nr             = epsilon_0 * m_e * (omega_r_smilei_SI**2) / (e**2) # reference density in SI
    Lr             = c / omega_r_smilei_SI                             # Smilei length unit in SI
    len_SI_fac     = lambda_r_SI / lamba_r_smilei                      # Length conversion factor from Smilei units to SI
    weight_SI_fac  = Nr * (Lr**3)                                      # Weight conversion factor from Smilei units to SI
    return time_SI_fac, len_SI_fac, weight_SI_fac

  def read_smilei_file(smilei_file_to_convert):
    # Read particle data from the sorted Smilei HDF5 file
    with h5py.File(smilei_file_to_convert, 'r') as file:
      x    = np.array(file['x'])
      y    = np.array(file['y'])
      z    = np.array(file['z'])
      px   = np.array(file['px'])
      py   = np.array(file['py'])
      pz   = np.array(file['pz'])
      w    = np.array(file['w'])
      time = np.array(file['Times'])
      time = time.astype(np.float64) # to interpret times as float not int
      return time, x, y, z, px, py, pz, w

  def get_number_of_particles(x):
    # Determine the total number of particles from Smilei simulation and display
    num_par = len(x[0,:])
    print('Total number of particles is: ' + str(num_par))
    return num_par

  def get_current_track(id, x, y, z, px, py, pz, w):
    # Extract current particle's trajectory
    x_curr    = x[:,id]
    y_curr    = y[:,id]
    z_curr    = z[:,id]
    px_curr   = px[:,id]
    py_curr   = py[:,id]
    pz_curr   = pz[:,id]
    w_curr    = w[:,id]
    return x_curr, y_curr, z_curr, px_curr, py_curr, pz_curr, w_curr

  def filter_out_nan(time, x_curr, y_curr, z_curr, px_curr, py_curr, pz_curr, w_curr):
    # Filter out NaN values in current particle's trajectory
    time_arr  = time[~np.isnan(x_curr)]
    x_arr     = x_curr[~np.isnan(x_curr)]
    y_arr     = y_curr[~np.isnan(y_curr)]
    z_arr     = z_curr[~np.isnan(z_curr)]
    px_arr    = px_curr[~np.isnan(px_curr)]
    py_arr    = py_curr[~np.isnan(py_curr)]
    pz_arr    = pz_curr[~np.isnan(pz_curr)]
    w_arr     = w_curr[~np.isnan(w_curr)]
    return time_arr, x_arr, y_arr, z_arr, px_arr, py_arr, pz_arr, w_arr

  def get_SI_unit_velocity(px_arr, py_arr, pz_arr):
    # Calculate SI unit velocity factor based on the simulation parameters
    gamma_arr = np.sqrt(1 + np.square(px_arr) + np.square(py_arr) + np.square(pz_arr)) 
    v_SI_fac  = c / gamma_arr 
    return v_SI_fac

  def print_writing_progress(id, par_fr):
    # Print progress every "par_fr" particles
    if float( id + 1) % par_fr == 0:
      print('Storing data from particle '+str( id + 1)+'.')      

  def create_new_file(time, x, y, z, px, py, pz, w,converted_file_name):
    # Create a new HDF5 file and store converted particle trajectory data
    with h5py.File(converted_file_name, 'w') as hdf_file:
      for id in range(0,num_par):

        x_curr, y_curr, z_curr, px_curr, py_curr, pz_curr, w_curr    = get_current_track(id, x, y, z, px, py, pz, w)
        time_arr, x_arr, y_arr, z_arr, px_arr, py_arr, pz_arr, w_arr = filter_out_nan(time, x_curr, y_curr, z_curr, px_curr, py_curr, pz_curr, w_curr)
        v_SI_fac                                                     = get_SI_unit_velocity(px_arr, py_arr, pz_arr)

        data_group = hdf_file.create_group(str(id+1))
        data_group.create_dataset('time', data = time_SI_fac * time_arr)
        data_group.create_dataset('x', data = len_SI_fac * x_arr )
        data_group.create_dataset('y', data = len_SI_fac * y_arr )
        data_group.create_dataset('z', data  = len_SI_fac * z_arr )    
        data_group.create_dataset('vx', data = v_SI_fac * px_arr)
        data_group.create_dataset('vy', data = v_SI_fac * py_arr)
        data_group.create_dataset('vz', data = v_SI_fac * pz_arr)  
        data_group.create_dataset('PIC_macroparticle_weight', data = weight_SI_fac * w_arr)  

        # Print progress every print_every particles
        print_writing_progress(id, print_every)

    print('File converted.')

  time_SI_fac, len_SI_fac, weight_SI_fac = get_conversion_factors(omega_r_smilei_SI, dt_smilei_SI)
  time, x, y, z, px, py, pz, w           = read_smilei_file(smilei_file_to_convert)
  num_par                                = get_number_of_particles(x)
  create_new_file(time, x, y, z, px, py, pz, w, converted_file_name)
