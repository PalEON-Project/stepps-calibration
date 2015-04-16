run1 = list(suff_fit = '12taxa_mid_comp_g_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = TRUE, EPs = FALSE)
run2 = list(suff_fit = '12taxa_mid_comp_vary_psi_EPs_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = TRUE, EPs = TRUE)
run3 = list(suff_fit = '12taxa_mid_comp_vary_psi_gamma_EPs_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = FALSE, EPs = TRUE)
run4 = list(suff_fit = '12taxa_mid_comp_vary_gamma_EPs_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = FALSE, EPs = TRUE)

runs = list(run1, run2, run3, run4)
runs=list(run3)
for (run in runs){
  source('r/cal_process.r')
}
