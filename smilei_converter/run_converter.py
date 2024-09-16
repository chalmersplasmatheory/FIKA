# Import required variables from the user input
from user_input_for_converter import *

# Import function that handles the conversion of simulation data from Smilei format to .h5 file suitable for the radiation calculation
from converter_from_smilei import convert_smilei_file

# This function call initiates the conversion process using the provided parameters and file names
convert_smilei_file(omega_r_smilei_SI, dt_smilei_SI, smilei_file_to_convert, converted_file_name, print_every)