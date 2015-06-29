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
There are several versions of the model implemented, each gets their own `.stan` file. Filenames with corresponding model description are below. Note thet `EPs` refers to the use of exchangeable priors.

  - `cal_g.stan` : base model, gaussian dispersal kernel
  - `cal_g_Kpsi_EPs.stan` : psi varies by taxon, gaussian dispersal kernel
  - `cal_g_Kgamma_EPs.stan` : gamma varies by taxon, gaussian dispersal kernel
  - `cal_g_Kpsi_Kgamma_EPs.stan` : psi and gamma vary by taxon, gaussian dispersal kernel
  - `cal_pl.stan` : base model, power-law kernel
  - `cal_pl_Ka_EPs.stan` :  a varies by taxon, power-law dispersal kernel
  - `cal_pl_Ka_Kb_EPs.stan` : a and b vary by taxon, power-law dispersal kernel
  - `cal_pl_Kgamma_EPs.stan` : gamma varies by taxon, power-law dispersal kernel
  - `cal_pl_Ka_Kgamma_EPs.stan` : a and gamma vary by taxon, power-law dispersal kernel
  - `cal_pl_Ka_Kb_Kgamma_EPs.stan` : a,b, and gamma vary by taxon, power-law dispersal kernel



  



