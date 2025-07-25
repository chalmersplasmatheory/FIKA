Script `get_trajectory_one_particle_transition regime.py` generates a HDF file with the trajectory of one electron with total (kinetic + rest) energy of $E_{e}=25$ MeV moving in the $x$ direction in plasma with plasma frequency $\omega_p$, where plasma density is $n_e=5\times 10^{18}$ $\rm{cm}^{-3}$. It oscillates in the $y$ direction with a simple sinusoidal motion with amplitude $r_\beta=0.5$ $\rm{\upmu m}$ and betatron frequency $\omega_\beta = \omega_p/\sqrt{2\gamma}$, where $\gamma$ is the electron Lorentz factor. The motion is in the transition regime between undulator and wiggler, i.e., undulator parameter $K = 1$. The on-axis fundamental frequency can be calculated as

$$
\omega_1 = \frac{2 \gamma^2 \omega_\beta}{1+K^2/2}, 
$$

where the undulator parameter corresponds to

$$
K = 1.33 \times 10^{-10} \sqrt{\gamma n_e [\rm{cm}^{-3}] } r_{\beta} [\rm{\upmu m}].
$$

The on-axis spectrum is expected to have peaks in the fundamental frequency $\omega_1\approx 26$ eV and several odd harmonics. After FIKA calculation, final data can be plotted with `plot_spectrum_one_particle_transition_regime.py`. This example can be compared with the results in Horný, V., et al. <em>Physics of Plasmas</em> 24.6 (2017): 063107.
#Note that the spectrum summation of one particle to generate file `final_spectrum.h5` was not necessary here, as the information was already obtained in `individual_spectra.h5`, but it was generated for testing purposes.

Files `final_spectrum.h5` and `individual_spectra.h5` were not included in the repository.

