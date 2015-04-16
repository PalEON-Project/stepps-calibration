// Author:  Andria Dawson
// Date:    October 2014
// Settlement era pollen estimation model based on STEPPS1
// Uses veg proportions and pollen counts to estimate process parameters
// With a Gaussian dispersal model using taxon-specific dispersal distance parameters psi[k]
 

data {
  int<lower=0> K;                // number of taxa
  int<lower=0> N_cores;          // number of core sites
  int<lower=0> N_cells;          // number of spatial cells
  int<lower=0> N_pot;            // number of potential contributing spatial cells

  int y[N_cores,K];              // pollen count data
  
  int idx_cores[N_cores];        // core cell indices
  
  vector[K] r[N_cells];          // composition proportions
  matrix[N_cells,N_cores] d;     // distance matrix
  matrix[N_pot, 2] d_pot;        // distances and counts of potential neighborhood
}
transformed data {
  matrix[N_cells,N_cores] d2;    // distance matrix squared
    
  d2 <- d .* d ;
}
parameters {
  vector<lower=0.01, upper=300>[K] phi;  // dirichlet precision pars
  real<lower=0.1, upper=2>    psi;       // dispersal par
  vector<lower=0, upper=1>[K] gamma;     // localness pars
  real<lower=-2, upper=2> mu_gamma;      // gamma hyperparameter
  real<lower=0> sigma_gamma;             // gamma hyperparameter
}

transformed parameters {
}

model {

  // declarations
  matrix[N_cells,N_cores] w;      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  real sum_w_pot;
  real max_r_new;
  int  max_r_new_idx;

  real N;
  real A;
  vector[K] alpha;
    
  // priors 
  phi         ~ uniform(0.01,300);
  psi         ~ uniform(0.1, 2);
  mu_gamma    ~ uniform(-2, 2);
  sigma_gamma ~ cauchy(0, 5);
 
  for (k in 1:K){
    increment_log_prob(- log(sigma_gamma) - log(gamma[k]) - log(1 - gamma[k]));
    increment_log_prob(- square(logit(gamma[k]) - mu_gamma) / (2 * square(sigma_gamma)));
  }  
  
  sum_w_pot <- 0;
  for (v in 1:N_pot)
    sum_w_pot <- sum_w_pot + d_pot[v,2] * exp(-square(d_pot[v,1])/square(psi));
  
  w <- exp(-(d2)/square(psi));
  
  for (i in 1:N_cores){
    for (k in 1:K){
      out_sum[k] <- 0;
      for (j in 1:N_cells){ // change N_hood to N_cells
	if (j != idx_cores[i]){
	  out_sum[k] <- out_sum[k] + w[j,i] * r[j][k];
	}  
      }
    }
    
    //local plus non-local piece
    for (k in 1:K)
      r_new[i,k] <-  gamma[k] * r[idx_cores[i],k] + out_sum[k] * (1-gamma[k]) / sum_w_pot;
    
    // // hacky!
    // // find taxon with highest proportional value
    // max_r_new <- 0;
    // for (k in 1:K){
    //   if (r_new[i,k] > max_r_new){
    //     max_r_new     <- r_new[i,k];
    //     max_r_new_idx <- k;
    //    }
    // }
    
    // for (k in 1:K){
    //   if (r_new[i,k] == 0){
    //     r_new[i,k] <- 0.0001;
    //     r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
    //     print("warning: zero proportion; core: ", i, "; taxon: ", k, " -> adjusting");
    //   }
    // }
  
    alpha <- phi .* r_new[i];
      
    A <- sum(alpha);
    N <- sum(y[i]);     

    increment_log_prob(lgamma(N + 1) + lgamma(A) - lgamma(N + A));
    for (k in 1:K) increment_log_prob( - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]));
    
  }
}
generated quantities{
  vector[N_cores] log_lik;

  {
    // declarations
    matrix[N_cells,N_cores] w;      
    vector[K] r_new[N_cores];
    vector[K] out_sum;    
  
    real sum_w_pot;

    real N;
    real A;
    vector[K] alpha;
    
    sum_w_pot <- 0;
    for (v in 1:N_pot)
      sum_w_pot <- sum_w_pot + d_pot[v,2] * exp(-square(d_pot[v,1])/square(psi));
  
    w <- exp(-(d2)/square(psi));
  
    for (i in 1:N_cores){
      for (k in 1:K){
	out_sum[k] <- 0;
	for (j in 1:N_cells){ // change N_hood to N_cells
	  if (j != idx_cores[i]){
	    out_sum[k] <- out_sum[k] + w[j,i] * r[j][k];
	  }  
	}
      }
    
      //local plus non-local piece
      for (k in 1:K)
	r_new[i,k] <-  gamma[k] * r[idx_cores[i],k] + out_sum[k] * (1-gamma[k]) / sum_w_pot;
      alpha <- phi .* r_new[i];
      
      A <- sum(alpha);
      N <- sum(y[i]);     

      log_lik[i] <- lgamma(N + 1) + lgamma(A) - lgamma(N + A);
      for (k in 1:K) log_lik[i] <- log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]);
    } 
  }
}