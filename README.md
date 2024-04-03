# FIKA
This repository contains the Far-field Intensity through Kinetic Analysis (FIKA) code, a Python-based computational tool designed for analyzing far-field radiation intensity from accelerated particles. FIKA leverages the Liénard-Wiechert potentials to calculate the radiation emitted by charged particles directly from their trajectories. The code calculates temporal profile of the radiated spectrum. It also employs Fourier transform techniques to analyze the frequency spectrum of the radiation.

To run the serial (nonparallel) version of FIKA, ensure you have the following Python packages installed, compatible with your version of Python:
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
    - `input_file` : Set the full path to the `.h5` file containing the information about the particle trajectories. Example: 
    ```bash
    input_file  = '/home/where_my_file_is/particle_trajectories.h5'
    ```
    The HDF5 file    contains groups represented as natural ascending numbers (in the string format), i. e., `'1'` , `'2'` , `'3'`  ,..... Each group represents one particle and contains datasets for the following quantities that need to be expressed in the SI units: 
      - `time` : Time is the units of s, a data of time values in which the particle trajectory was tracked. 
      - `x` : x position of the particle in the units of m, expressed in corresponding times
      - `y` : y position of the particle in the units of m, expressed in corresponding times
      - `z` : z position of the particle in the units of m, expressed in corresponding times
      - `vx` : velocity of the particle in the x direction in the units of m/s, expressed in corresponding times
      - `vy` : velocity of the particle in the y direction in the units of m/s, expressed in corresponding times
      - `vz` : velocity of the particle in the z direction in the units of m/s, expressed in corresponding times
      - `PIC_macroparticle_weight` (OPTIONAL) : If the particle data are obtained from a particle-in-cell (PIC) code, one "particle" typically corresponds to a macroparticle, representing a real number of electrons. PIC_macroparticle_weight is this real number of particles represented by the macroparticle. It needs to be included if the calculated radiation unit is required to be realistic.  
        
    - `output_folder`: Path where the results will be saved. Example:

   
   ```bash
    output_folder  = '/home/where_my_output_will_be/'
    ```
      
      
    - `PIC_macroparticle_weights`: Takes boolean `True` or `False`. Use `True` only if you included `PIC_macroparticle_weight` parameter in the input file and you want to include the macropaticle weights (real number of electrons representing the macroparticle) in the calculations.
      
    - `E_radMax_eV` : Maximum radiation energy expected in the units of eV. The sampling of Fourier transformation will be 20x this value. The mmaximum energy has to be higher than 2x this value for correct sampling. If you find out that the highest photon energy is higher than this value, we recommend to rerun the code increasing this value. Example:
    
    ```bash
    E_radMax_eV   = 1000 
    ```
    
    - `E_slice_eV` : A bin energy size in the units of eV of the final output energy spectrum histogram. Example:

    
    ```bash     
    E_slice_eV    = 2 
    ```
    
    - `t_slice_s` : A bin time size in the units of second of the final output of the radiation temporal profile. Example:

    
    ```bash     
    t_slice_s     = 5e-18
    ```
    
    - `phi` : Angles φ in the spherical coordinates, giving the position of the observer. Example:
    
    ```bash     
    phi           = 0
    ```
    
    - `theta` : Angle θ in the spherical coordinates, giving the position of the observer. Example:
    
    ```bash     
    theta         = pi/2
    ```
    
4. Run the simulation using the following command:

    ```bash
    python run.py
    ```
<!--
## Authors
FIKA is authored by Dominika Maslarova and the members of the
[Chalmers Plasma Theory group](https://ft.nephy.chalmers.se/).
-->
