from input import *
from calculateRadiation import calculate_spectrum_one_particle, calculate_spectrum_all_particles
import os 

# Get the name of the variable
def get_variable_name(variable):
    for name in globals():
        if id(variable) == id(globals()[name]):
            return name
    return None

# Check if the variables are either an integer or a float
def check_if_the_value_numeric(variable):
  if not isinstance(variable, (int, float)):
    print("Value of", variable, "is not numeric.")
    # Stop the program
    raise SystemExit('Stopping the program because variable '+ str(get_variable_name(variable))+ ' is not numeric.')

# Check if the value is boolean
def check_if_the_value_boolean(variable):
  if not isinstance(variable, bool):
    print("Value of", variable, "should be either True or False.")
    # Stop the program
    raise SystemExit('Stopping the program because variable '+ str(get_variable_name(variable))+ ' is not boolean.')

# Check if the file exists
def check_if_the_file_exists(file_path):
  if not os.path.exists(file_path):
    print("The file path",file_path,"does not exist.")
    # Stop the program
    raise SystemExit("Stopping the program because the input file was not found.")

# Check if the folder exists
def check_if_the_folder_exists(folder_path):
  if not os.path.isdir(folder_path):
    print("The output path",folder_path,"does not exist.")
    # Stop the program
    raise SystemExit("Stopping the program because the output directory was not found.")

# Execute this part when running the script: 
if __name__ == "__main__":
  # Check the user input
  check_if_the_value_numeric(charge)
  check_if_the_value_numeric(E_radMax_eV)
  check_if_the_value_numeric(r)
  check_if_the_value_numeric(phi)
  check_if_the_value_numeric(theta)
  check_if_the_value_boolean(PIC_macroparticle_weights)
  check_if_the_file_exists(input_file)
  check_if_the_folder_exists(output_folder)
  check_if_the_value_boolean(sum_spectra)
  check_if_the_value_numeric(print_every_par_spectrum)
  if sum_spectra == True:
    check_if_the_value_numeric(E_slice_eV)
    check_if_the_value_numeric(t_slice_s)
    check_if_the_value_numeric(print_every_spectrum_sum)

  # Run the simulation  
  if parallel == 'Disabled':
    calculate_spectrum_one_particle(charge, E_radMax_eV, r, phi, theta, input_file, output_folder, PIC_macroparticle_weights, print_every_par_spectrum)
    if sum_spectra == True:
      calculate_spectrum_all_particles(output_folder, E_slice_eV, t_slice_s, PIC_macroparticle_weights, print_every_spectrum_sum)
    else:
      print('Calculation finished.')
  else:
    print('Parallel version is not yet available. Choose parallel = \'Disabled\' for serial calculation.')

