Script `get_trajectory_one_particle_undulator.py` generates a HDF file with the trajectory of one electron with total (kinetic + rest) energy of $E_{e}=10$ MeV moving in the $x$ direction in plasma with plasma frequency $\omega_p$, where plasma density is $n_e=5\times 10^{24}$ $\rm{cm}^{-3}$. It oscilates in the $y$ direction with a simple sinusoidal motion with amplitude $r_\beta=0.05$ $\rm{\upmu m}$ and betatron frequency $\omega_\beta = \omega_p/\sqrt{2\gamma}$, where $\gamma$ is the electron Lorentz factor. The motion is undulator-like, i.e. undulator parameter $K \ll 1$. The axis spectrum is expected to be almost monoenergetic with fundamental frequency $\omega_1$. After FIKA calculation, final data can be plotted with `plot_spectrum_one_particle_undulator.py`. The result can be compared with analytical formula 

$$
\omega_1 = \frac{2 \gamma^2 \omega_\beta}{1+K^2/2} = 10.2~\rm{eV},
$$

where the undulator parameter can be calculated as

$$
K = 1.33 \times 10^{-10} \sqrt{\gamma n_e} r_{\beta}.
$$

Similar example can be found in the dissertation of Vojtěch Horný, CTU in Prague, Czech Republic, 2018.
