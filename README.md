# stepps-calibration
==================

This repository contains the STEPPS pollen-veg calibration model code and paper.

To run the code:
  1.  Install Stan and RStan.
  2.  Edit the pointer to Stan in the Makefile.
  3.  Run make. To compile `calibration.stan`, type: `make calibration.exe` in the terminal.
  4.  Run the code. I'll make a script to do this so it's easy (now that our data set is fixed).

## model versions
==================
There are several versions of the model implemented, each gets their own `.stan` file. Filenames with corresponding model description are below.

  - `calibration.stan` : base model with gaussian dispersal kernel
  - `calibration_pl.stan` : base model with power-law kernel
  - `calibration_vary_gamma.stan` : gamma varies by taxon, gaussian dispersal kernel
  - `calibration_vary_gamma_psi.stan` : gamma and psi both vary by taxon
  - `calibration_vary_gamma_EPs.stan` : gamma varies by taxon, with exchangeable priors
  - `calibration_vary_gamma_psi_EPs` : gamma and psi both vary by taxon, with exchangeable priors


  



