# compute effective sample size form a stanfit object
ess <- function(fit){
  ess = summary(fit)$summary[,"n_eff"]
  return(ess)
}

# compute acceptance rate from a stanfit object
ar <- function(fit){
  post = extract(fit, permuted=FALSE, inc_warmup=FALSE)
  ar = apply(post, "chains", FUN = function(x) nrow(unique(as.data.frame(x)))) / nrow(post) # acceptance rates
  return(ar)
}

get_phi_stats <- function(phi, taxa){
  phi_bar <- colMeans(phi)
  phi_L   <- apply(phi, 2, function(x){quantile(x, probs = (0.025))})
  phi_U   <- apply(phi, 2, function(x){quantile(x, probs = (0.975))})
  
  phi_stats <- data.frame(mean = phi_bar, L = phi_L, U = phi_U, taxa = taxa)
  
  phi_stats
}


compute_C <- function(post, N_pot, d_pot, kernel){
  
  col_substr = substr(colnames(post[,1,]), 1, 3)
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = mean(post[,1,which(col_substr == 'psi')])
    } else {
      psi   = colMeans(post[,1,which(col_substr == 'psi')])
    }
    
    C = sum( d_pot[,2] * exp(-d_pot[,1]^2/psi^2) )
  } else if (kernel=='pl'){
    a   = mean(post[,1,which(col_substr == 'a')])
    b   = mean(post[,1,which(col_substr == 'b')])
    C = sum( d_pot[,2] * (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d_pot[,1] / a) )
  }
  
  
  return(C)
  
}

#predicted pollen based on weighted neighborhoods using estimated pars
pollen_preds <- function(post, N_cores, d, idx_cores, r, C, one_psi, kernel){
  
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_substr = substr(colnames(post[,1,]), 1, 3)
  
  #   phi   = summary(fit)$summary[,'mean'][1:K]
  #   psi   = summary(fit)$summary[,'mean'][K+1]
  #   gamma = summary(fit)$summary[,'mean'][K+2]
  
  phi   = colMeans(post[,1,which(col_substr == 'phi')])
  gamma = mean(post[,1,which(col_substr == 'gam')])
  
  #   
  sum_w = vector(length=N_cores, mode='numeric')
  r_new = matrix(NA, nrow=N_cores, ncol=K)
  preds = matrix(NA, nrow=N_cores, ncol=K)
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = mean(post[,1,which(col_substr == 'psi')])
    } else {
      psi   = colMeans(post[,1,which(col_substr == 'psi')])
    }
    if (one_psi){
      w = exp(-(d*d)/(psi*psi))
    } else{
      w = array(NA, c(K, N_cells, N_cores))
      for (k in 1:K){ 
        w[k,,] = exp(-(d*d)/(psi[k]*psi[k]))
      }
    }
  } else if (kernel=='pl'){
    print("Kernel type : (inverse) power law")
    a = mean(post[,1,which(col_substr == 'a')])
    b = mean(post[,1,which(col_substr == 'b')])
    w = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a) 
  }
  
  for (i in 1:N_cores){
    print(i)
    
    #       for (j in 1:N_hood){
    #         if (idx_hood[i,j] > 0){
    #           w[i,j] <- exp(-(d[idx_hood[i,j], i])^2/(psi)^2)
    #         } 
    #       }
    
    out_sum = rep(0, K)
    #     sum_w <- 0
    
    if (one_psi){
      for (j in 1:N_cells){ # changed N_hood to N_locs
        if (j != idx_cores[i]){
          #         if (idx_hood[i,j] > 0){
          #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
          out_sum <- out_sum + w[j,i]*r[j,]
        }  
      }
    } else {
      for (k in 1:K){
        print(paste0('k = ', k))
        for (j in 1:N_cells){ # changed N_hood to N_locs
          if (j != idx_cores[i]){
            #         if (idx_hood[i,j] > 0){
            #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
            out_sum[k] <- out_sum[k] + w[k,j,i]*r[j,k]
          }  
        }
      }
      
    }
    #     sum_w   <- sum(out_sum)
    #     print(sum_w)
    
    sum_w[i] = sum(out_sum)
    #     r_new[i,]  = gamma*r[idx_cores[i],] + (1-gamma)*out_sum/sum_w
    r_new[i,]  = gamma*r[idx_cores[i],] + (1-gamma)*out_sum/C
    preds[i,] = phi*r_new[i,]    
    
  }
  
  alpha = rowSums(preds)   
  
  #convert to proportions
  preds = t(apply(preds, 1, function(x) x/sum(x)))
  
  return(list(preds=preds, alpha=alpha, sum_w=sum_w))
}

#predicted pollen based on weighted neighborhoods using estimated pars
pollen_preds_distance <- function(post, N_cores, d, idx_cores, r, C, radius){
  
  rescale = 1e6
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_substr = substr(colnames(post[,1,]), 1, 3)
  
  phi   = colMeans(post[,1,which(col_substr == 'phi')])
  gamma = mean(post[,1,which(col_substr == 'gam')])
  if (one_psi){
    psi   = mean(post[,1,which(col_substr == 'psi')])
  } else {
    psi   = colMeans(post[,1,which(col_substr == 'psi')])
  }
  
  r_local = matrix(NA, nrow=N_cores, ncol=K)
  r_nl = matrix(NA, nrow=N_cores, ncol=K)
  r_int = matrix(NA, nrow=N_cores, ncol=K)
  r_new = matrix(NA, nrow=N_cores, ncol=K)
  preds = matrix(NA, nrow=N_cores, ncol=K)
  preds_tot = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_int = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_dist = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_dist_taxon = array(NA, c(K, N_cores, length(radius)))
  #   w     = matrix(0, nrow=N_cores, ncol=N_cells)
  
  w = exp(-(d*d)/(psi*psi))
  
  for (rad in 1:length(radius)){
    print(paste0('rad = ', rad))
    for (i in 1:N_cores){
      idx_int = which(d[,i]<radius[rad]/rescale)
      
      #       for (j in 1:N_hood){
      #         if (idx_hood[i,j] > 0){
      #           w[i,j] <- exp(-(d[idx_hood[i,j], i])^2/(psi)^2)
      #         } 
      #       }
      
      out_sum = rep(0, K)
      out_sum_int = rep(0, K)
      sum_w <- 0
      
      for (j in 1:N_cells){ # changed N_hood to N_locs
        if (j != idx_cores[i]){
          #         if (idx_hood[i,j] > 0){
          #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
          out_sum <- out_sum + w[j,i]
          
          if (j %in% idx_int){
            out_sum_int <- out_sum_int + w[j,i]
          }
        }  
      }
      
      sum_w   <- sum(out_sum)
      
      sum_w_int   <- sum(out_sum_int)
      
      #       r_local[i,] = r[idx_cores[i],]
      r_nl[i,]  = gamma + (1-gamma) / C * out_sum
      r_int[i,] = gamma + (1-gamma) / C * out_sum_int
      #       r_new[i,]  = gamma*r[idx_cores[i],] + (1-gamma)*out_sum/sum_w
      #       preds[i,] = phi*r_new[i,]  
      
    }
    
    # do i want to multiply by phi??
    #     preds_loc = rowSums(r_local)
    preds_tot[rad,]  = r_nl[,1]#rowSums(r_nl)
    preds_int[rad,] = r_int[,1]#rowSums(r_int)
    preds_dist[rad,] = preds_int[rad,]/preds_tot[rad,] 
    preds_dist_taxon[,,rad] = r_int/r_nl 
    
  }
  
  #   convert to proportions
  #   preds = t(apply(preds, 1, function(x) x/sum(x)))
  
  return(list(preds_tot=preds_tot, preds_int=preds_int, preds_dist=preds_dist))
}

#predicted pollen based on weighted neighborhoods using estimated pars
dispersal_decay <- function(post, d_pot, C, radius, kernel){
  
  rescale = 1e6
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_substr = substr(colnames(post[,1,]), 1, 3)
  
  phi   = colMeans(post[,1,which(col_substr == 'phi')])
  gamma = mean(post[,1,which(col_substr == 'gam')])
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = mean(post[,1,which(col_substr == 'psi')])
    } else {
      psi   = colMeans(post[,1,which(col_substr == 'psi')])
    }
    w      = exp(-(d_pot*d_pot)/(psi*psi))
  } else if (kernel=='pl'){
    a = mean(post[,1,which(col_substr == 'a')])
    b = mean(post[,1,which(col_substr == 'b')])
    w = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d_pot / a) 
  }
  
  r_int  = vector(length=length(radius), mode='numeric')
  C_test = sum(w)
  
  for (rad in 1:length(radius)){
    idx_int    = which(d_pot[,1]<radius[rad]/rescale)
    sum_w_int  = sum(w[idx_int,1])
    r_int[rad] = gamma + (1-gamma) / C * sum_w_int     
  }
  
  return(r_int)
}


#scale veg props in focal cell based on phi
phi_scale_veg <- function(post, N_cores, r, idx_cores){
  
  iters = dim(post)[1]
  K     = dim(r)[2]
  
  phi   = summary(fit)$summary[,'mean'][1:K]
  
  r_new   = matrix(NA, nrow=N_cores, ncol=K)
  phi_veg = matrix(NA, nrow=N_cores, ncol=K)
  
  for (i in 1:N_cores){
    for (k in 1:K){
      
      r_new[i,k] = phi[k]*r[idx_cores[i],k]
      
    }
  }
  
  preds = r_new/apply(r_new, 1, sum) 
  
  return(preds)
  
}

get_quants <- function(fit, npars){
  quants <- cbind(summary(fit)$summary[,'mean'][1:npars],
                  summary(fit)$summary[,'2.5%'][1:(npars)],
                  summary(fit)$summary[,'50%'][1:(npars)],
                  summary(fit)$summary[,'97.5%'][1:(npars)],
                  summary(fit)$summary[,'n_eff'][1:(npars)],
                  summary(fit)$summary[,'Rhat'][1:(npars)])
  colnames(quants) = c('mean', '2.5%', '50%', '97.5%', 'n_eff', 'Rhat')
  
  return(quants)
}

compute_props <- function(y, taxa){
  pollen_props <- y/rowSums(y) 
  colnames(pollen_props) <- taxa
  
  return(pollen_props)
}

power_law <- function(d, a, b) {
  x = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a)
 return(x)
} 
