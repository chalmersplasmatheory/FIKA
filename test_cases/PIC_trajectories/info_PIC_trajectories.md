This is an example of electron trajectories obtained from a particle-in-cell simulation with the SMILEI code [1]. The file for the FIKA input `PIC_particle_set.h5` was obtained with the Smilei converter (see `README.md`). The radiation has synchrotron-like character, corresponding to betatron radiation from laser wakefield acceleration. Electron beam is moving in the positive $x$ direction undergoing small transverse betatron oscillations in the $y$ & $z$ directions. The `final_spectrum.h5` file contains the spectrum summation of the whole electron beam. The `individual_spectra.h5` file containing the individual spectra of particles is not shared here git but was produced during the calculation.

The macroparticle spectra (`PIC_macroparticle_weights = False` option in `input.py`) output can be found in the `without_macroparticle_weights/` folder.  The `final_spectrum.h5` file contains the spectrum summation of the whole electron beam. The `individual_spectra.h5` file containing the individual spectra of particles is not shared on git, but was produced during the calculation.

Real-particle spectrum (`PIC_macroparticle_weights = True` option in `input.py`) can be found in the `with_macroparticle_weights/final_spectrum.h5`.
Both cases can be plotted with `plot_spectrum_PIC_trajectories.py` for comparison.

[1] J. Derouillat, A. Beck, F. Pérez, T. Vinci, M. Chiaramello, A. Grassi, M. Flé, G. Bouchard, I. Plotnikov, N. Aunai, J. Dargent, C. Riconda and M. Grech, SMILEI: a collaborative, open-source, multi-purpose particle-in-cell code for plasma simulation, Comput. Phys. Commun. 222, 351-373 (2018), arXiv:1702.05128

