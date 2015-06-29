g = list(suff_fit  = 'cal_g_ALL_v0.3', 
             suff_dat = '12taxa_mid_comp_ALL_v0.2',
             kernel    = 'gaussian', 
             one_psi   = TRUE, 
             one_gamma = TRUE, 
             EPs       = FALSE)
g_Kpsi = list(suff_fit  = 'cal_g_Kpsi_EPs_ALL_v0.3', 
              suff_dat = '12taxa_mid_comp_ALL_v0.2',
              kernel    = 'gaussian', 
              one_psi   = FALSE, 
              one_gamma = TRUE, 
              EPs       = TRUE)
g_Kpsi_Kgamma = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_ALL_v0.3',
                     suff_dat = '12taxa_mid_comp_ALL_v0.2',
                     kernel    = 'gaussian', 
                     one_psi   = FALSE, 
                     one_gamma = FALSE, 
                     EPs       = TRUE)
g_Kgamma = list(suff_fit  = 'cal_g_Kgamma_EPs_ALL_v0.3', 
                suff_dat = '12taxa_mid_comp_ALL_v0.2',
                kernel    = 'gaussian', 
                one_psi   = TRUE, 
                one_gamma = FALSE, 
                EPs       = TRUE)
pl = list(suff_fit  = 'cal_pl_ALL_v0.3', 
                  suff_dat = '12taxa_mid_comp_ALL_v0.2',
                  kernel    = 'pl', 
                  one_a     = TRUE,
                  one_b     = TRUE,
                  one_gamma = TRUE, 
                  EPs       = FALSE)
pl_Kgamma = list(suff_fit  = 'cal_pl_Kgamma_EPs_ALL_v0.3', 
                 suff_dat = '12taxa_mid_comp_ALL_v0.2',
                 kernel    = 'pl', 
                 one_a     = TRUE,
                 one_b     = TRUE,
                 one_gamma = FALSE, 
                 EPs       = TRUE)
pl_Ka = list(suff_fit  = 'cal_pl_Ka_EPs_ALL_v0.3', 
                 suff_dat = '12taxa_mid_comp_ALL_v0.2',
                 kernel    = 'pl', 
                 one_a     = FALSE,
                 one_b     = TRUE,
                 one_gamma = TRUE, 
                 EPs       = TRUE)
pl_Ka_Kb = list(suff_fit  = 'cal_pl_Ka_Kb_EPs_ALL_v0.3', 
             suff_dat = '12taxa_mid_comp_ALL_v0.2',
             kernel    = 'pl', 
             one_a     = FALSE,
             one_b     = FALSE,
             one_gamma = TRUE, 
             EPs       = TRUE)
pl_Ka_Kb_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kb_Kgamma_EPs_ALL_v0.3', 
                       suff_dat = '12taxa_mid_comp_ALL_v0.2',
                       kernel    = 'pl', 
                       one_a     = FALSE,
                       one_b     = FALSE,
                       one_gamma = FALSE, 
                       EPs       = TRUE)
pl_Ka_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3', 
                suff_dat = '12taxa_mid_comp_ALL_v0.2',
                kernel    = 'pl', 
                one_a     = FALSE,
                one_b     = TRUE,
                one_gamma = FALSE, 
                EPs       = TRUE)


# runs = list(g, g_Kpsi, g_Kgamma, g_Kpsi_Kgamma)
runs = list(pl, pl_Ka, pl_Ka_Kb, pl_Ka_Kb_Kgamma, pl_Kgamma, pl_Ka_Kgamma)
runs = list(pl, pl_Ka, pl_Ka_Kb, pl_Ka_Kb_Kgamma, pl_Kgamma, pl_Ka_Kgamma)
runs = list(pl_Ka_Kb, pl_Ka_Kb_Kgamma, pl_Kgamma, pl_Ka_Kgamma)
runs = list(pl_Ka_Kgamma)
for (run in runs){
  source('r/cal_process.r')
}
