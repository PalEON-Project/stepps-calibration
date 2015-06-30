g = list(suff_fit  = 'cal_g_ALL_v0.4', 
         kernel    = 'gaussian', 
         one_psi   = TRUE, 
         one_gamma = TRUE, 
         EPs       = FALSE,
         handle    = 'g')
pl = list(suff_fit  = 'cal_pl_ALL_v0.4', 
          kernel    = 'pl', 
          one_a     = TRUE, 
          one_b     = TRUE,
          one_gamma = TRUE, 
          EPs       = FALSE,
          handle    = 'pl')

g_Kpsi_Kgamma = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_ALL_v0.4', 
                     kernel    = 'gaussian', 
                     one_psi   = FALSE, 
                     one_gamma = FALSE, 
                     EPs       = TRUE,
                     handle    = 'g_Kpsi_Kgamma')
pl_Ka_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.4', 
                    kernel    = 'pl', 
                    one_a     = FALSE, 
                    one_b     = TRUE,
                    one_gamma = FALSE, 
                    EPs       = TRUE,
                    handle    = 'pl_Ka_Kgamma')
runs = list(g, pl, g_Kpsi_Kgamma, pl_Ka_Kgamma)