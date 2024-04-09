# FIKA

This repository contains the Far-field Intensity through Kinetic Analysis (FIKA) code, a Python-based computational tool designed for analyzing far-field radiation intensity from accelerated particles. FIKA leverages the Liénard-Wiechert potentials to calculate the radiation emitted by charged particles directly from their trajectories. The code calculates the temporal profile of the radiated spectrum. It also uses Fourier transform techniques to analyze the frequency spectrum of the radiation.

## Requirements

Python 3 is required to run the FIKA code, with version 3.7 or higher recommended for optimal compatibility and performance. To run the serial (non-parallel) version of FIKA, ensure you have the following Python packages installed, compatible with your version of Python:
- `h5py`
- `matplotlib`
- `numpy`
- `scipy`

An MPI-supported parallel version is under development and will be available soon.

## Running FIKA

To execute the FIKA simulation, follow these steps:

1. Navigate to the directory containing the FIKA code.

2. Modify the parameters in `input.py` according to your simulation requirements.
  The following parameters need to be specified in the code:
    - `input_file`
          Set the full path to the `.h5` file containing the information about the particle trajectories. Example:
      
          input_file = '/home/where_my_file_is/particle_trajectories.h5'
      
        The HDF5 file contains groups represented as natural ascending numbers (in the string format), i. e., `'1'` , `'2'` , `'3'`  ,..... Each group represents one particle and contains datasets for the following quantities that need to be expressed in the SI units: 
       - `time` : A dataset of time values (in seconds) at which the particle trajectory was tracked.
       - `x` : The x position of the particle (in meters), corresponding to each time value.
       - `y` : The y position of the particle (in meters), corresponding to each time value.
       - `z` : The z position of the particle (in meters), corresponding to each time value.
       - `vx` : The velocity of the particle in the x direction (in meters per second), corresponding to each time value.
       - `vy` : The velocity of the particle in the y direction (in meters per second), corresponding to each time value.
       - `vz` : The velocity of the particle in the z direction (in meters per second), corresponding to each time value.
       - `PIC_macroparticle_weight` (OPTIONAL) : If the particle data are obtained from a particle-in-cell (PIC) code, one "particle" typically corresponds to a macroparticle, representing a real number of particles. PIC_macroparticle_weight is this real number of particles represented by the macroparticle. It needs to be included if the calculated radiation unit is required to be realistic.  

   

    - `output_folder`
    
         Path where the results will be saved. Example:
    
           output_folder = '/home/where_my_output_will_be/'
           
    
    - `PIC_macroparticle_weights`
    
       Accepts a boolean value, `True` or `False`. Set to `True` if you have included the PIC_macroparticle_weight parameter in the input file and wish to incorporate the macroparticle weights (the real number of particles represented by each macroparticle) in the calculations. Example:
    
           PIC_macroparticle_weights = True

    - `parallel`
    
       The only option currently available is 'Disabled'. Support for parallel options is under development. Example:    
       
           parallel = 'Disabled'
      
    - `charge`
    
       Specifies the charge of one particle. In the case of PIC macroparticles, specify the charge of the actual particle represented. Example (for an electron particle/macroparticle):      
            from scipy.constants import e
            charge = -e  
      
      
    - `E_radMax_eV`
    
         Defines the maximum radiation energy expected, in units of eV. The Fourier transformation sampling will be 20 times this value. The actual maximum energy must be higher than twice this value for accurate sampling. If the user discovers that the highest photon energy exceeds this value, it is recommended to rerun the code with an increased `E_radMax_eV`. Example: 
    
          E_radMax_eV = 1000 
   
    - `r`

       Specifies the radial position in spherical coordinates (in meters), indicating the observer's position. Example:
      
          r = 1


    - `phi`

       Specifies azimuthal angle φ in spherical coordinates (in radians), indicating the observer's position. Example:
      
          phi = 0

    
    - `theta`
      
        Specifies polar angle θ in spherical coordinates (in radians), indicating the observer's position. Example:

          from numpy import pi
          theta = pi/2

    - `print_every_par_spectrum`

       Determines how often the information about the calculation progress will be written to the standard output. The specified number indicates the number of new particles processed after which the output will be printed. Example:

          print_every_par_spectrum  = 200 

    -  `sum_spectra`
    
        If `True`, the energy spectra and temporal profiles of individual particles will be summed to provide compact information about the entire particle beam. This summation process begins after the radiation from all particles has been calculated. If `False`, the calculation is finished right after the individual spectra are calculated. Example:

          sum_spectra  = True


        The following additional parameters need to be specified in the output file with `sum_spectra  = True` option: 
            
       - `E_slice_eV` 
      
            Specifies the bin energy size in eV for the final output energy spectrum histogram. Example:

            E_slice_eV = 2 

       - `t_slice_s`
      
            Specifies the bin time size in seconds for the final output of the radiation temporal profile. Example:
   
            t_slice_s = 5e-18

       -  `print_every_spectrum_sum` 
       
            Indicates how often information about the summation progress will be written to the standard output. The number specifies how many new particles have been summed at each update. Example:

            print_every_spectrum_sum  = 200                                        


4. Run the simulation using the following command:

    ```bash
    python main.py
    ```

## Output 

The output in the form of HDF file, named `individual_spectra.h5`, saved in the folder specified by the user. The HDF5 file contains groups represented as natural ascending numbers (in the string format), i. e., `'1'` , `'2'` , `'3'`, ... representing the individual particles' IDs corresponding to the original `input_file`. Each group contains following datasets: 
   - `ene` : Energies of radiated photons in eV.
   - `spectrum_ene` : Photon energy distribution $\\frac{d^2 W}{d E d\Omega}\$ received by the observer, expressed in the units of J/eV/sr.
   - `t` : Range of the observer time in seconds.
   - `spectrum_t` : Spectral intensity $\\frac{d^2 I}{d t d\Omega}\$ corresponding to the observer time in the units of J/s/sr.

In case `sum_spectra = True` ia specified by the user, another HDF5 file, `final_spectrum.h5`, is created, cointaining following groups:
   - `ene` : Energies of radiated photons in eV.
   - `spectrum_ene` : Photon energy distribution \(\frac{d^2 W}{d E d\Omega}\) received by the observer summed for all the particles, expressed in the units of J/eV/sr.
   - `t` : Range of the observer time in s.
   - `spectrum_t` : Spectral intensity \(\frac{d^2 I}{d t d\Omega}\) summed for all the particles, corresponding to the observer time in the units of J/s/sr.


## Smilei converter

For users of the Smilei particle-in-cell code [1], an extra module which converts the Smilei output into a FIKA input file is available in the `smilei_converter` folder. To run the converter, follow these steps:

1. Navigate to the directory containing the FIKA code.

2. Navigate to `smilei_converter/`

3. Modify the parameters in `user_input_for_converter.py`. The following parameters need to be specified:

    - `omega_r_smilei_SI/`

        Specifies reference angular frequency from your Smilei simulation in the SI units. When running the Smilei simulations, all the Smilei units are derived from the reference angular frequency. Be aware that the parameter here requires the value of the reference frequency in the SI units! In the following example, the frequency is set to corresponding to 0.8 micron wavelength:

        omega_r_smilei_SI      = 2354564459136066.5

    - `dt_smilei_SI/`

        Specifies the PIC timestep used in your Smilei simulation (i.e. the timestep of the PIC loop). It is gain necessary to convert this value to the SI units! In the following example, timestep is equal to 8.539240837072692e-17 s:

        dt_smilei_SI           = 8.539240837072692e-17       

    - `smilei_file_to_convert/`

        Indicates the full path to the Smilei .h5 file containing particle tracking data, including the file name. It has to be the file that is already sorted by postprocessing with happi library. This can be done in TrackParticles diagnostic by choosing `sort = True` in TrackParticles(..., sort=True, ...). See the Smilei documentation [2], "Open a TrackParticles diagnostic", for more information. Unsorted version cannot be processed! The file needs to contain the following information about the particles need to be tracked: `'x', 'y', 'z', 'px', 'py', 'pz', 'w'`. Example for a sorted track file with species "electron":

        smilei_file_to_convert = 'TrackParticles_electron.h5' 

    - `converted_file_name/`

        Indicates the full path for the output file .h5 where the data for radiation calculation in FIKA will be stored. Example:
         
        converted_file_name    = 'test_particle_set_1.h5'   

    - `converted_file_name/`
        
        Progress update frequency. After a new set of "print_every" particles is processed, the update will be printed into a standard output. Example, where the update is printed after the data of each 100 new particles are converted:

        print_every            = 100

### References
[1]  J. Derouillat, A. Beck, F. Pérez, T. Vinci, M. Chiaramello, A. Grassi, M. Flé, G. Bouchard, I. Plotnikov, N. Aunai, J. Dargent, C. Riconda and M. Grech, SMILEI: a collaborative, open-source, multi-purpose particle-in-cell code for plasma simulation, Comput. Phys. Commun. 222, 351-373 (2018), arXiv:1702.05128

[2]  [Smilei documentation on post-processing](https://smileipic.github.io/Smilei/Use/post-processing.html) Accessed April 9, 2024.

<!--
## Authors
FIKA is authored by Dominika Maslarova and the members of the
[Chalmers Plasma Theory group](https://ft.nephy.chalmers.se/).
-->
