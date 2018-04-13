// Author:  Andria Dawson
// Date:    October 2014
// Settlement era pollen estimation model based on STEPPS1
// Uses veg proportions and pollen counts to estimate process parameters
// With an Inverse Power Law dispersal kernel with taxon-specific a and gamma
 

data {
  int<lower=0> K;                // number of taxa
  int<lower=0> N_cores;          // number of core sites
  int<lower=0> N_cells;          // number of spatial cells
  int<lower=0> N_pot;            // number of potential contributing spatial cells

  int y[N_cores,K];              // pollen count data
  
  int idx_cores[N_cores];        // core cell indices
  int N_hood[N_cores];            // hood cell counts by core
  int idx_hood[N_cores, N_cells]; // hood cell indices
  
  vector[K] r[N_cells];          // composition proportions
  matrix[N_cells,N_cores] d;     // distance matrix 
  matrix[N_pot, 2] d_pot;        // distances and counts of potential neighborhood

  vector[K] mu_a;
  vector[K] sigma_a;
}
transformed data {
 vector[N_cores] N;
  vector[N_cores] lgamma_Nplus1;

 for (i in 1:N_cores){
   N[i] <- sum(y[i]);
   lgamma_Nplus1[i] <- lgamma(N[i] + 1);
 }
}
parameters {
  vector<lower=0.01, upper=300>[K] phi;  // dirichlet precision pars
  vector<lower=0, upper=1>[K] gamma;     // localness pars

  real<lower=-2, upper=2> mu_gamma;  // gamma hyperparameter
  real<lower=0> sigma_gamma;             // gamma hyperparameter

  vector<lower=log(1e-4), upper=log(1)>[K] log_a;
  real<lower=2.001, upper=6> b;

  //real<lower=log(1e-4), upper=log(1)> mu_a;
  //real<lower=1e-5> sigma_a;
}
transformed parameters {
  vector[K] a;

  for (k in 1:K){
    a[k] <- exp(log_a[k]);
  }
}

model {

  // declarations
  matrix[N_cells,N_cores] w[K];      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  vector[K] sum_w_pot;
  vector[K] kernel_p1;
  vector[K] log_kernel_p1;

  // real max_r_new;
  // int  max_r_new_idx;
    
  real A;
  vector[K] alpha;

  // priors 
  phi     ~ uniform(0.01,300);
  mu_gamma    ~ uniform(-2, 2);
  sigma_gamma ~ cauchy(0, 5); 
  b ~ uniform(2.001, 6);  

  for (k in 1:K){
    log_a[k] ~ normal(mu_a[k], sigma_a[k]);
    increment_log_prob(- log(sigma_gamma) - log(gamma[k]) - log(1 - gamma[k]));
    increment_log_prob(- square(logit(gamma[k]) - mu_gamma) / (2 * square(sigma_gamma)));
  }  
  
  for (k in 1:K){
    //kernel_p1[k] <-  (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] );
    log_kernel_p1[k] <-  log(b-2) + log(b-1) - log(2 * pi()) - 2*log(a[k]);
    sum_w_pot[k] <- 0; 
    for (v in 1:N_pot){
      //sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * kernel_p1[k] * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
      sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * exp(log_kernel_p1[k]  - b *  log(1 + d_pot[v,1] / a[k]) ); 
    }
    for (i in 1:N_cells)
      for (j in 1:N_cores)
	// w[k][i,j] <-  kernel_p1[k] * pow( 1 + d[i,j] / a[k], -b[k]) ;
	w[k][i,j] <-  log_kernel_p1[k] - b * log( 1 + d[i,j] / a[k]) ;

    w[k] <- exp(w[k]); 
  }
  
  for (i in 1:N_cores){        
    for (k in 1:K){
      out_sum[k] <- 0;
      for (j in 1:N_hood[i]){
	out_sum[k] <- out_sum[k] + w[k][idx_hood[i,j],i]*r[idx_hood[i,j]][k];
      }  
    }

    //local plus non-local piece
    for (k in 1:K)
      r_new[i,k] <- gamma[k]*r[idx_cores[i],k] + out_sum[k] * (1-gamma[k])/sum_w_pot[k];
        
    // // when zeros in raw data, readjust to non-zero
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
    //}
      
    alpha <- phi .* r_new[i];
    A     <- sum(alpha);     

    increment_log_prob(lgamma_Nplus1[i] + lgamma(A) - lgamma(N[i] + A));
    for (k in 1:K) increment_log_prob( - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]));
  }
}
generated quantities {
  vector[N_cores] log_lik;

  {
    // declarations
    matrix[N_cells,N_cores] w[K];      
    vector[K] r_new[N_cores];
    vector[K] out_sum;    
  
    vector[K] sum_w_pot;
    vector[K] kernel_p1;
    vector[K] log_kernel_p1;

    // real max_r_new;
    // int  max_r_new_idx;
    
    real A;
    vector[K] alpha;
  
    for (k in 1:K){
      //kernel_p1[k] <-  (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] );
      log_kernel_p1[k] <-  log(b-2) + log(b-1) - log(2 * pi()) - 2*log(a[k]);
      sum_w_pot[k] <- 0; 
      for (v in 1:N_pot){
	//sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * kernel_p1[k] * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
	sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * exp(log_kernel_p1[k]  - b *  log(1 + d_pot[v,1] / a[k]) ); 
      }
      for (i in 1:N_cells)
	for (j in 1:N_cores)
	  // w[k][i,j] <-  kernel_p1[k] * pow( 1 + d[i,j] / a[k], -b[k]) ;
	  w[k][i,j] <-  log_kernel_p1[k] - b * log( 1 + d[i,j] / a[k]) ;

      w[k] <- exp(w[k]); 
    }
  
    for (i in 1:N_cores){        
      for (k in 1:K){
	out_sum[k] <- 0;
	// for (j in 1:N_cells){ // change N_hood to N_cells
	//   if (j != idx_cores[i]){
	//     out_sum[k] <- out_sum[k] + w[k][j,i]*r[j][k];
	//   }  
	// }
	for (j in 1:N_hood[i]){
	  out_sum[k] <- out_sum[k] + w[k][idx_hood[i,j],i]*r[idx_hood[i,j]][k];
	}  
      }

      //local plus non-local piece
      for (k in 1:K)
	r_new[i,k] <- gamma[k]*r[idx_cores[i],k] + out_sum[k]*(1-gamma[k])/sum_w_pot[k];
      
      alpha <- phi .* r_new[i];
      A     <- sum(alpha);     
    
      log_lik[i] <- lgamma_Nplus1[i] + lgamma(A) - lgamma(N[i] + A);
      for (k in 1:K) log_lik[i] <- log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]);
    
    }
  } 
}
