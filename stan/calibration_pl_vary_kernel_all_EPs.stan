// Author:  Andria Dawson
// Date:    October 2014
// Settlement era pollen estimation model based on STEPPS1
// Uses veg proportions and pollen counts to estimate process parameters
// With a Gaussian dispersal model using taxon-specific dispersal distance parameters psi[k]
 

data {
  int<lower=0> K;                // number of taxa
  int<lower=0> N_cores;          // number of core sites
  int<lower=0> N_cells;          // number of spatial cells
  int<lower=0> N_pot;           // number of potential contributing spatial cells


  int y[N_cores,K];              // pollen count data
  
  int idx_cores[N_cores];        // core cell indices
  
  vector[K] r[N_cells];          // pls proportions
  matrix[N_cells,N_cores] d;    // distance matrix squared
  matrix[N_pot, 2] d_pot;        // distances and counts of potential neighborhood
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
  // real<lower=0, upper=1> gamma;          // local proportion par
 vector<lower=0, upper=1>[K] gamma;

  real<lower=-100, upper=100> mu_gamma;       // psi hyperparameter
  real<lower=0> sigma_gamma;             // psi hyperparameter

  vector<lower=log(1e-6), upper=log(500)>[K] log_a;
  vector<lower=log(2), upper=log(100)>[K] log_b;

  real<lower=log(1e-6), upper=log(500)> mu_a;
  real<lower=1e-6> sigma_a;
  real<lower=log(2), upper=log(100)> mu_b;
  real<lower=1e-6> sigma_b;
}

transformed parameters {
  vector[K] a;
  vector[K] b;

  for (k in 1:K){
    a[k] <- exp(log_a[k]);
    b[k] <- exp(log_b[k]);
  }

}

model {

  // declarations
  matrix[N_cells,N_cores] w[K];      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  //real sum_w;
  vector[K] sum_w_pot;
  vector[K] kernel_p1;
  vector[K] log_kernel_p1;

  real max_r_new;
  int  max_r_new_idx;

  // priors 
  phi     ~ uniform(0.01,300);
  mu_gamma    ~ uniform(-10, 10);
  sigma_gamma ~ cauchy(0, 10);
  // gamma   ~ uniform(0,1); 
  mu_a    ~ uniform(log(1e-6), log(500));
  sigma_a ~ cauchy(1e-6, 2);  
  mu_b    ~ uniform(log(2), log(100));
  sigma_b ~ cauchy(1e-6, 2);  
  
  for (k in 1:K){
    log_a[k] ~ normal(mu_a, sigma_a);
    log_b[k] ~ normal(mu_b, sigma_b);

    increment_log_prob(- log(sigma_gamma) - log(gamma[k]) - log(1 - gamma[k]));
    increment_log_prob(- square(logit(gamma[k]) - mu_gamma) / (2 * square(sigma_gamma)));
  }  
  
  for (k in 1:K){
    //kernel_p1[k] <-  (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] );
    log_kernel_p1[k] <-  log(b[k]-2) + log(b[k]-1) - log(2 * pi()) - 2*log(a[k]);

    sum_w_pot[k] <- 0;
  
    for (v in 1:N_pot){
      //sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * kernel_p1[k] * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
      sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * exp(log_kernel_p1[k]  - b[k] *  log(1 + d_pot[v,1] / a[k]) ); 
    }
    for (i in 1:N_cells)
      for (j in 1:N_cores)
	// w[k][i,j] <-  kernel_p1[k] * pow( 1 + d[i,j] / a[k], -b[k]) ;
	w[k][i,j] <-  log_kernel_p1[k] - b[k] * log( 1 + d[i,j] / a[k]) ;

    w[k] <- exp(w[k]);
      
  }
  
  for (i in 1:N_cores){
    // for (j in 1:N_hood){
    //   if (idx_hood[i,j] > 0){
    //     w[i,j] <- exp(-square(d[i,idx_hood[i,j]]/psi));
    //   } else {
    //     w[i,j] <- 0;
    //   }
    //  }
   
   
   
    // if (sum(r[idx_cores[i]]) == 0 ){
    //   print("Cell ", idx_cores[i], " has NO vegetation! Neighborhood props will not sum to 1.");
    // }
   
    
    // local piece
    for (k in 1:K)
      r_new[i] <- gamma[k]*r[idx_cores[i]];
   
    for (k in 1:K) {out_sum[k] <- 0;}
    //sum_w <- 0;
    
    for (k in 1:K){
      for (j in 1:N_cells){ // change N_hood to N_cells
	if (j != idx_cores[i]){
	  out_sum[k] <- out_sum[k] + w[k][j,i]*r[j][k];
	}  
      }
    }

    //sum_w   <- sum(out_sum);
     
    //print("out_sum", sum(out_sum));
    //print("sum_w", sum_w);
    
    //local vs. non-local
    for (k in 1:K)
      r_new[i,k] <- r_new[i,k] + out_sum[k]*(1-gamma[k])/sum_w_pot[k];
    
    //print(r_new[i]);
    //print(sum(r_new[i]));
    
    /*
    if (sum(r_new[i]) < 1-1e-6 ){
      print("i ", i);
      print("core props ", r[idx_cores[i]]);
      print("sum r_new ", sum(r_new[i]));
    }
    */
    
    // hacky!
    // find taxon with highest proportional value
    max_r_new <- 0;
    for (k in 1:K){
      if (r_new[i,k] > max_r_new){
        max_r_new     <- r_new[i,k];
        max_r_new_idx <- k;
       }
    }
    
    for (k in 1:K){
      if (r_new[i,k] == 0){
        r_new[i,k] <- 0.0001;
        r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
        print("warning: zero proportion; core: ", i, "; taxon: ", k, " -> adjusting");
      }
    }
    
    
    {
      real A;
      vector[K] alpha;
      
      alpha <- phi .* r_new[i];
      A     <- sum(alpha);     

      print("A = ", A);

      increment_log_prob(lgamma_Nplus1[i] + lgamma(A) - lgamma(N[i] + A));
      for (k in 1:K) increment_log_prob( - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]));
    }
    
  }
}
generated quantities {
  vector[N_cores] log_lik;

  {
  // declarations
  matrix[N_cells,N_cores] w[K];      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  //real sum_w;
  vector[K] sum_w_pot;
  vector[K] kernel_p1;

  real max_r_new;
  int  max_r_new_idx;
 
  for (k in 1:K){
    kernel_p1[k] <-  (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] );

    sum_w_pot[k] <- 0;
  
    for (v in 1:N_pot)
      sum_w_pot[k] <- sum_w_pot[k] + d_pot[v,2] * kernel_p1[k] * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
 
    for (i in 1:N_cells)
      for (j in 1:N_cores)
	w[k][i,j] <- kernel_p1[k] * pow( 1 + d[i,j] / a[k], -b[k]) ;
  }
  
  for (i in 1:N_cores){
    // for (j in 1:N_hood){
    //   if (idx_hood[i,j] > 0){
    //     w[i,j] <- exp(-square(d[i,idx_hood[i,j]]/psi));
    //   } else {
    //     w[i,j] <- 0;
    //   }
    //  }
   
   
   
    // if (sum(r[idx_cores[i]]) == 0 ){
    //   print("Cell ", idx_cores[i], " has NO vegetation! Neighborhood props will not sum to 1.");
    // }
   
    
    // local piece
    for (k in 1:K)
      r_new[i] <- gamma[k]*r[idx_cores[i]];
   
    for (k in 1:K) {out_sum[k] <- 0;}
    //sum_w <- 0;
    
    for (k in 1:K){
      for (j in 1:N_cells){ // change N_hood to N_cells
	if (j != idx_cores[i]){
	  out_sum[k] <- out_sum[k] + w[k][j,i]*r[j][k];
	}  
      }
    }

    //sum_w   <- sum(out_sum);
     
    //print("out_sum", sum(out_sum));
    //print("sum_w", sum_w);
    
    //local vs. non-local
    for (k in 1:K)
      r_new[i,k] <- r_new[i,k] + out_sum[k]*(1-gamma[k])/sum_w_pot[k];
    
    //print(r_new[i]);
    //print(sum(r_new[i]));
    
    /*
    if (sum(r_new[i]) < 1-1e-6 ){
      print("i ", i);
      print("core props ", r[idx_cores[i]]);
      print("sum r_new ", sum(r_new[i]));
    }
    */
    
    // hacky!
    // find taxon with highest proportional value
    max_r_new <- 0;
    for (k in 1:K){
      if (r_new[i,k] > max_r_new){
        max_r_new     <- r_new[i,k];
        max_r_new_idx <- k;
       }
    }
    
    for (k in 1:K){
      if (r_new[i,k] == 0){
        r_new[i,k] <- 0.0001;
        r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
        print("warning: zero proportion; core: ", i, "; taxon: ", k, " -> adjusting");
      }
    }
    
    
    {
      real A;
      vector[K] alpha;
      
      alpha <- phi .* r_new[i];
      A     <- sum(alpha);     

      //print("A = ", A);

      log_lik[i] <- lgamma_Nplus1[i] + lgamma(A) - lgamma(N[i] + A);
      for (k in 1:K) log_lik[i] <- log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]);
    }
    
  }
  } 
}