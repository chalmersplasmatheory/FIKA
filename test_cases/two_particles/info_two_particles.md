Script `get_trajectory_two_particles_transition regime.py` generates 2 HDF files `final_spectrum.h5` into folders `one/` and `two/`. 

The file in folder `one/` corresponds to the trajectory of one electron with total (kinetic + rest) energy of $E_{e}=150$ MeV moving in the $x$ direction in plasma with plasma frequency $\omega_p$, where plasma density is $n_e=5\times 10^{18}$ $\rm{cm}^{-3}$. It oscillates in the $y$ direction with a simple sinusoidal motion with amplitude $r_\beta=1.2$ $\rm{\upmu m}$ and betatron frequency $\omega_\beta = \omega_p/\sqrt{2\gamma}$, where $\gamma$ is the electron Lorentz factor. This motion corresponds to the wiggler regime, which typically occurs in betatron radiation during laser wakefield acceleration.
The same example of radiation spectrum can be found in Horný, V., et al. <em>Physics of Plasmas</em> 24.6 (2017): 063107.

File `final_spectrum.h5` in folder `two/`corresponds to the spectrum resulting from the summation of two identical electron trajectories, as described above. This artificial test case is designed to verify whether the spectrum summation operates correctly; summing two identical particles should result in doubled intensity values.

After completing the FIKA calculation, the final data can be visualized using `plot_spectrum_two_particles.py`. 