# FIKA

This repository contains the Far-field Intensity through Kinetic Analysis (FIKA) code, a Python-based computational tool designed for analyzing far-field radiation intensity from accelerated particles. FIKA leverages the Liénard-Wiechert potentials to calculate the radiation emitted by charged particles directly from their trajectories. The code calculates temporal profile of the radiated spectrum. It also employs Fourier transform techniques to analyze the frequency spectrum of the radiation.

## Requirements

Python 3 is required to run the FIKA code (version 3.7 and higher recommended). To run the serial (nonparallel) version of FIKA, ensure you have the following Python packages installed, compatible with your version of Python:
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
       - `time` : Time is the units of s, a data of time values in which the particle trajectory was tracked. 
       - `x` : x position of the particle in the units of m, expressed in corresponding times
       - `y` : y position of the particle in the units of m, expressed in corresponding times
       - `z` : z position of the particle in the units of m, expressed in corresponding times
       - `vx` : velocity of the particle in the x direction in the units of m/s, expressed in corresponding times
       - `vy` : velocity of the particle in the y direction in the units of m/s, expressed in corresponding times
       - `vz` : velocity of the particle in the z direction in the units of m/s, expressed in corresponding times
       - `PIC_macroparticle_weight` (OPTIONAL) : If the particle data are obtained from a particle-in-cell (PIC) code, one "particle" typically corresponds to a macroparticle, representing a real number of electrons. PIC_macroparticle_weight is this real number of particles represented by the macroparticle. It needs to be included if the calculated radiation unit is required to be realistic.  

            


    - `output_folder`
    
         Path where the results will be saved. Example:
    
           output_folder = '/home/where_my_output_will_be/'
           
    
    - `PIC_macroparticle_weights`
    
       Takes boolean `True` or `False`. Use `True` only if you included `PIC_macroparticle_weight` parameter in the input file and you want to include the macropaticle weights (real number of electrons representing the macroparticle) in the calculations. Example:
    
           PIC_macroparticle_weights = True

    - `parallel`
    
       Only option currently available is `'Disabled'`. Parallel options coming soon. Example:
    
           parallel = 'Disabled'
      
    - `charge`
    
       Charge of one particle. In case of macroparticles, also put the charge of the REAL particle. Example (electron particle/macroparticle):
      
            from scipy.constants import e
            charge = -e  
      
      
    - `E_radMax_eV`
      
         Maximum radiation energy expected in the units of eV. The sampling of Fourier transformation will be 20x this value. The maximum energy has to be higher than 2x this value for correct sampling. If you find out that the highest photon energy is higher than this value, we recommend to rerun the code increasing this value. Example:
    
    
          E_radMax_eV = 1000 
   
    - `r`

       Radial position in the spherical coordinates in the units of m, giving the position of the observer. Example:
      
          r = 1


    - `phi`

       Angle φ in the spherical coordinates, giving the position of the observer. Example:
      
          phi = 0

    
    - `theta`
      
        Angle θ in the spherical coordinates, giving the position of the observer. Example:

          from numpy import pi
          theta = pi/2

    - `print_every_par_spectrum`

       Specifies how often the information about the calculation progress will be written into a standard output. The number specifies the new number of particles processed after which the output will be printed. Example:

          print_every_par_spectrum  = 200 

    -  `sum_spectra`
    
        If `True`, the energy spectra and temporal profiles of individual particles will be summed to obtain a compact information about the whole particle beam. The summation process starts after radiation from all the particles is calculated. If `False`, the calculation is finished right after the individual spectra are calculated.

        Following additional parameters need to be specified in the output file 
            
       - `E_slice_eV` 
      
        A bin energy size in the units of eV of the final output energy spectrum histogram. Example:

            E_slice_eV = 2 

       - `t_slice_s`
      
         A bin time size in the units of second of the final output of the radiation temporal profile. Example:
   
            t_slice_s = 5e-18

       -  `print_every_spectrum_sum` 
       
         Specifies how often the information about the summation progress will be written into a standard output. The number specifies the number of    particles that are already summed. Example:

            print_every_spectrum_sum  = 200                                        


4. Run the simulation using the following command:

    ```bash
    python run.py
    ```

## Output 

The output in the form of HDF file `individual_spectra.h5`, saved in the folder specified by the user. The HDF5 file contains groups represented as natural ascending numbers (in the string format), i. e., `'1'` , `'2'` , `'3'`, ... representing the individual particles IDs corresponding to the original `input_file`. Each group contains following datasets: 
   - `freq` : Energies of radiated photons in eV.
   - `spectrum_freq` : Corresponding spectral intensity for each photon energy $\(\frac{d^2 I}{d E d\Omega}\)$, expressed in the units of J/eV/sr.
   - `t` : Range of the observer time in s.
   - `spectrum_t` : Spectral intensity $\(\frac{d^2 I}{d t d\Omega}\)$ corresponding to the observer time in the units of J/s/sr.

In case `sum_spectra=True` ia specified by the user, another HDF5 file `final_spectrum.h5` is created, cointaining following groups:
   - `freq` : Energies of radiated photons in eV.
   - `spectrum_freq` : Corresponding spectral intensity for each photon energy \(\frac{d^2 I}{d E d\Omega}\) summed for all the particles, expressed in the units of J/eV/sr.
   - `t` : Range of the observer time in s.
   - `spectrum_t` : Spectral intensity \(\frac{d^2 I}{d t d\Omega}\) summed for all the particles, corresponding to the observer time in the units of J/s/sr.


<!--
## Authors
FIKA is authored by Dominika Maslarova and the members of the
[Chalmers Plasma Theory group](https://ft.nephy.chalmers.se/).
-->
