# run1 = list(suff_fit = '12taxa_mid_comp_g_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = TRUE, EPs = FALSE)
# run2 = list(suff_fit = '12taxa_mid_comp_vary_psi_EPs_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = TRUE, EPs = TRUE)
# run3 = list(suff_fit  = '12taxa_mid_comp_vary_psi_gamma_EPs_v0.3', 
#             kernel    = 'gaussian', 
#             one_psi   = FALSE, 
#             one_gamma = FALSE, 
#             EPs       = TRUE)
# run4 = list(suff_fit = '12taxa_mid_comp_vary_gamma_EPs_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = FALSE, EPs = TRUE)
# 
# run = list(suff_fit  = 'test_vary_psi_gamma_EPs_WAIC', 
#            kernel    = 'gaussian', 
#            one_psi   = FALSE, 
#            one_gamma = FALSE, 
#            EPs       = TRUE)

run1 = list(suff_fit  = 'cal_g_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = TRUE, 
            one_gamma = TRUE, 
            EPs       = FALSE)
run2 = list(suff_fit  = 'cal_g_Kpsi_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = TRUE, 
            EPs       = TRUE)
run3 = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = FALSE, 
            EPs       = TRUE)
run4 = list(suff_fit  = 'cal_g_Kgamma_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = TRUE, 
            one_gamma = FALSE, 
            EPs       = TRUE)
run5 = list(suff_fit  = 'cal_pl_v0.3', 
            kernel    = 'pl', 
            one_a     = TRUE,
            one_b     = TRUE,
            one_gamma = TRUE, 
            EPs       = FALSE)
run6 = list(suff_fit  = 'cal_pl_Kgamma_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = TRUE,
            one_b     = TRUE,
            one_gamma = FALSE, 
            EPs       = TRUE)
run7 = list(suff_fit  = 'cal_pl_Ka_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = FALSE,
            one_b     = TRUE,
            one_gamma = TRUE, 
            EPs       = TRUE)
run8 = list(suff_fit  = 'cal_pl_Ka_Kb_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = FALSE,
            one_b     = FALSE,
            one_gamma = TRUE, 
            EPs       = TRUE)
run9 = list(suff_fit  = 'cal_pl_Ka_Kb_Kgamma_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = FALSE,
            one_b     = FALSE,
            one_gamma = FALSE, 
            EPs       = TRUE)
run10 = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = FALSE,
            one_b     = TRUE,
            one_gamma = FALSE, 
            EPs       = TRUE)


run_g = list(suff_fit  = 'cal_g_UMW_v0.3', 
             suff_dat = '12taxa_mid_comp_UMW_v0.2',
             kernel    = 'gaussian', 
             one_psi   = TRUE, 
             one_gamma = TRUE, 
             EPs       = FALSE)
run_pl = list(suff_fit  = 'cal_pl_UMW_v0.3',
              suff_dat = '12taxa_mid_comp_UMW_v0.2',
              kernel    = 'pl', 
              one_a     = TRUE,
              one_b     = TRUE,
              one_gamma = TRUE, 
              EPs       = FALSE)

g_all = list(suff_fit  = 'cal_g_ALL_v0.3', 
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
pl_all = list(suff_fit  = 'cal_pl_ALL_v0.3', 
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

runs = list(run1, run2, run3, run4, run5, run6)
runs = list(run3, run4, run5, run6, run7, run8, run9)

# runs = list(run_g, run_pl, run_g_all, run_pl_all)

runs = list(g_Kpsi, g_Kgamma, g_Kpsi_Kgamma, pl_Kgamma)
for (run in runs){
  source('r/cal_process.r')
}
