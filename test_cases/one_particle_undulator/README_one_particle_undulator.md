Script `get_trajectory_one_particle_undulator.py` generates a HDF file with the trajectory of one electron with total (kinetic + rest) energy of $E_{e}=10$ MeV moving in the $x$ direction in plasma with plasma frequency $\omega_p$, where plasma density is $n_e=5\times 10^{18}$ $\rm{cm}^{-3}$. It oscillates in the $y$ direction with a simple sinusoidal motion with amplitude $r_\beta=0.05$ $\rm{\upmu m}$ and betatron frequency $\omega_\beta = \omega_p/\sqrt{2\gamma}$, where $\gamma$ is the electron Lorentz factor. The motion is undulator-like, i.e. undulator parameter $K \ll 1$. The axis spectrum is expected to be almost monoenergetic with fundamental frequency $\omega_1 \approx$ 10 eV. After FIKA calculation, final data can be plotted with `plot_spectrum_one_particle_undulator.py`. The result can be compared with analytical formula 

$$
\omega_1 = \frac{2 \gamma^2 \omega_\beta}{1+K^2/2},
$$

where the undulator parameter can be calculated as

$$
K = 1.33 \times 10^{-10} \sqrt{\gamma n_e [\rm{cm}^{-3}]} r_{\beta} [\rm{\upmu m}].
$$

Similar example can be found in the dissertation of Vojtěch Horný, CTU in Prague, Czech Republic, 2018.

Note that the spectrum summation of one particle to generate file `final_spectrum.h5` was not necessary here, as the information was already obtained in `individual_spectra.h5`, but it was generated for testing purposes.

Files `final_spectrum.h5` and `individual_spectra.h5` were not included in the repository.
