# 
# #predicted pollen based on weighted neighborhoods using estimated pars
# pollen_preds <- function(post, N_cores, d, idx_cores, r, sum_w, one_psi, one_gamma, kernel){
#   
#   iters   = dim(post)[1]
#   K       = ncol(r)
#   N_cells = nrow(d)
#   
#   col_substr = substr(colnames(post[,1,]), 1, 3)
#   
#   #   phi   = summary(fit)$summary[,'mean'][1:K]
#   #   psi   = summary(fit)$summary[,'mean'][K+1]
#   #   gamma = summary(fit)$summary[,'mean'][K+2]
#   #   
#   phi   = colMeans(post[,1,which(col_substr == 'phi')])
#   #   gamma = mean(post[,1,which(col_substr == 'gam')])
#   
#   #   
#   #   sum_w = vector(length=N_cores, mode='numeric')
#   r_new = matrix(NA, nrow=N_cores, ncol=K)
#   preds = matrix(NA, nrow=N_cores, ncol=K)
#   log_lik = array(NA, c(N_cores))
#   
#   if (one_gamma){
#     gamma = rep(mean(post[,1,which(col_substr == 'gam')]), K)
#   } else {
#     gamma = colMeans(post[,1,which(col_substr == 'gam')])
#   }
#   
#   if (kernel=='gaussian'){
#     if (one_psi){
#       psi   = rep(mean(post[,1,which(col_substr == 'psi')]), K)
#     } else {
#       psi   = colMeans(post[,1,which(col_substr == 'psi')])
#     }
#     #     w = array(NA, c(K, N_cells, N_cores))
#     #     for (k in 1:K){ 
#     #       w[k,,] = exp(-(d*d)/(psi[k]*psi[k]))
#     #     }
#     
#   } else if (kernel=='pl'){
#     print("Kernel type : (inverse) power law")
#     a = mean(post[,1,which(col_substr == 'a')])
#     b = mean(post[,1,which(col_substr == 'b')])
#     w = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a) ^ (-b)
#   }
#   
#   for (i in 1:N_cores){
#     print(i)
#     
#     #       for (j in 1:N_hood){
#     #         if (idx_hood[i,j] > 0){
#     #           w[i,j] <- exp(-(d[idx_hood[i,j], i])^2/(psi)^2)
#     #         } 
#     #       }
#     
#     out_sum = rep(0, K)
#     #     sum_w <- 0
#     
#     #     if (one_psi){
#     #       for (j in 1:N_cells){ # changed N_hood to N_locs
#     #         if (j != idx_cores[i]){
#     #           #         if (idx_hood[i,j] > 0){
#     #           #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
#     #           out_sum <- out_sum + w[j,i]*r[j,]
#     #         }  
#     #       }
#     #     } else {
#     for (k in 1:K){
#       print(paste0('k = ', k))
#       for (j in 1:N_cells){ # changed N_hood to N_locs
#         if (j != idx_cores[i]){
#           #         if (idx_hood[i,j] > 0){
#           #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
#           #     for (k in 1:K){ 
#           #       w[k,,] = exp(-(d*d)/(psi[k]*psi[k]))
#           #     }
#           if (kernel == 'gaussian'){
#             w = exp(-(d[j,i]*d[j,i])/(psi[k]*psi[k]))
#             out_sum[k] <- out_sum[k] + w*r[j,k]
#           } else if (kernel == 'pl'){
#             out_sum[k] <- out_sum[k] + w[j,i]*r[j,k]
#             #             out_sum[k] <- out_sum[k] + w[k,j,i]*r[j,k]
#           }
#         }  
#       }
#     }
#     
#     #     }
#     #     sum_w   <- sum(out_sum)
#     #     print(sum_w)
#     
#     #     sum_w[i] = sum(out_sum)
#     #     r_new[i,]  = gamma*r[idx_cores[i],] + (1-gamma)*out_sum/sum_w
#     
#     alpha = rep(NA, K)
#     for (k in 1:K){
#       r_new[i,k]  = gamma[k]*r[idx_cores[i],k] + (1-gamma[k])*out_sum[k]/sum_w[k]
#       alpha[k] = phi[k]*r_new[i,k]        
#     }
#     
#     A = sum(alpha)   
#     N = sum(y[i,])
#     
#     log_lik[i] = lgamma(N + 1) + lgamma(A) - lgamma(N + A)
#     for (k in 1:K)  log_lik[i] = log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]) 
#     
#   }
#   
# 
#   
#   return(list(log_lik))
# }


#predicted pollen based on weighted neighborhoods using estimated pars
log_lik_iter <- function(post, N_cores, d, idx_cores, r, sum_w, one_psi, one_gamma, kernel){
  
  iter = 1
  
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_substr = substr(colnames(post[,1,]), 1, 3)

  phi   = post[iter,1,which(col_substr == 'phi')]
 
  r_new = matrix(NA, nrow=N_cores, ncol=K)
  preds = matrix(NA, nrow=N_cores, ncol=K)
  log_lik = array(NA, c(N_cores))
  
  if (one_gamma){
    gamma = rep(post[iter,1,which(col_substr == 'gam')], K)
  } else {
    gamma = post[iter,1,which(col_substr == 'gam')]
  }
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = rep(post[iter,1,which(col_substr == 'psi')], K)
    } else {
      psi   = post[iter,1,which(col_substr == 'psi')]
    }
    
  } else if (kernel=='pl'){
    print("Kernel type : (inverse) power law")
    a = mean(post[,1,which(col_substr == 'a')])
    b = mean(post[,1,which(col_substr == 'b')])
    w = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a) ^ (-b)
  }
  
  sum_w=rep(0, K)
  for (k in 1:K)
    sum_w[k] = sum(d_pot[,2] * exp(-d_pot[,1]^2/psi[k]^2))
  
  for (i in 1:N_cores){
    print(i)

    out_sum = rep(0, K)
    for (k in 1:K){
      print(paste0('k = ', k))
      for (j in 1:N_cells){ # changed N_hood to N_locs
        if (j != idx_cores[i]){
          if (kernel == 'gaussian'){
            w = exp(-(d[j,i]*d[j,i])/(psi[k]*psi[k]))
            out_sum[k] <- out_sum[k] + w*r[j,k]
          } else if (kernel == 'pl'){
            out_sum[k] <- out_sum[k] + w[j,i]*r[j,k]
          }
        }  
      }
    }
    
    alpha = rep(NA, K)
    for (k in 1:K){
      r_new[i,k]  = gamma[k]*r[idx_cores[i],k] + (1-gamma[k])*out_sum[k]/sum_w[k]
    }
    
    max_r_new_idx = which.max(r_new[i,])
    
    for (k in 1:K){
      if (r_new[i,k] == 0){
        r_new[i,k] <- 0.0001;
        r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
        print(paste0("warning: zero proportion; core: ", i, "; taxon: ", k, " -> adjusting"))
      }
    }
    
    alpha = phi*r_new[i,]        
    A = sum(alpha)   
    N = sum(y[i,])
    
    log_lik[i] = lgamma(N + 1) + lgamma(A) - lgamma(N + A)
    for (k in 1:K)  log_lik[i] = log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]) 
    
  }
  
  return(list(log_lik))
}