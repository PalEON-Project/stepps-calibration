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
  matrix[N_cells,N_cores] d2;    // distance matrix squared
  matrix[N_pot, 2] d_pot;        // distances and counts of potential neighborhood
}

transformed data {
}

parameters {
  vector<lower=0.01, upper=300>[K] phi;  // dirichlet precision pars
  real<lower=0.1, upper=2>    psi;  // dispersal par
  //vector<lower=0, upper=1>[K] gamma;     // local proportion par

  //vector<lower=log(0.1), upper=log(2)>[K] log_psi;
  //vector<lower=log(1e-15), upper=log(1)>[K] log_gamma;
  vector<lower=0, upper=1>[K] gamma;

  //real<lower=log(0.1), upper=log(2)> mu_psi;         // psi hyperparameter
  //real<lower=0> sigma_psi;               // psi hyperparameter

  real<lower=-100, upper=100> mu_gamma;       // psi hyperparameter
  real<lower=0> sigma_gamma;             // psi hyperparameter
}

transformed parameters {
  //  vector[K] psi;
  //vector[K] gamma;

  // for (k in 1:K){
  //   psi[k]   <- exp(log_psi[k]); 
    //gamma[k] <- exp(log_gamma[k]); 
  // }
}

model {

  // declarations
  matrix[N_cells,N_cores] w;      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  real sum_w;
  real sum_w_pot;
  real max_r_new;
  int  max_r_new_idx;

  // priors 
  phi         ~ uniform(0.01,300);
  psi         ~ uniform(0.1, 2);
  mu_gamma    ~ uniform(-10, 10);
  sigma_gamma ~ cauchy(0, 10);
  // mu_psi      ~ uniform(log(0.1), log(2));
  // sigma_psi   ~ cauchy(0, 4);
  //phi         ~ uniform(0.01,300);
  //gamma     ~ uniform(0,1); 
  for (k in 1:K){
    //psi[k] ~ uniform(0.1,2);
    //log_gamma[k] ~ normal(mu_gamma, sigma_gamma);
    increment_log_prob(- log(sigma_gamma) - log(gamma[k]) - log(1 - gamma[k]));
    increment_log_prob(- square(logit(gamma[k]) - mu_gamma) / (2 * square(sigma_gamma)));
    //log_psi[k]   ~ normal(mu_psi, sigma_psi);
  }  
  
  sum_w_pot <- 0;
  for (v in 1:N_pot)
    sum_w_pot <- sum_w_pot + d_pot[v,2] * exp(-square(d_pot[v,1])/square(psi));
  //print("C = ", C);  
  
  w <- exp(-(d2)/square(psi));
  
  for (i in 1:N_cores){
 /*   for (j in 1:N_hood){
      if (idx_hood[i,j] > 0){
        w[i,j] <- exp(-square(d[i,idx_hood[i,j]]/psi));
      } else {
        w[i,j] <- 0;
      }
     }
   */ 
   
   /*
    if (sum(r[idx_cores[i]]) == 0 ){
      print("Cell ", idx_cores[i], " has NO vegetation! Neighborhood props will not sum to 1.");
    }
    */
    
    // local piece
    for (k in 1:K)
      r_new[i] <- gamma[k] * r[idx_cores[i]];
   
    for (k in 1:K) {out_sum[k] <- 0;}
    sum_w <- 0;
    
    for (k in 1:K){
      for (j in 1:N_cells){ // change N_hood to N_cells
	if (j != idx_cores[i]){
	  out_sum[k] <- out_sum[k] + w[j,i] * r[j][k];
	}  
      }
    }

    //sum_w   <- sum(out_sum);
     
    //print("out_sum", sum(out_sum));
    //print("sum_w", sum_w);
    
    //local vs. non-local
    for (k in 1:K)
      r_new[i,k] <- r_new[i,k] + out_sum[k] * (1-gamma[k]) / sum_w_pot;
    
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
      real N;
      real A;
      vector[K] alpha;
      
      alpha <- phi .* r_new[i];
      
      A <- sum(alpha);
      N <- sum(y[i]);     

      increment_log_prob(lgamma(N + 1) + lgamma(A) - lgamma(N + A));
      for (k in 1:K) increment_log_prob( - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]));
    }
    
  }
}
generated quantities{
  vector[N_cores] log_lik;

  {
    // declarations
    matrix[N_cells,N_cores] w;      
    vector[K] r_new[N_cores];
    vector[K] out_sum;    
  
    real sum_w;
    real sum_w_pot;
    real max_r_new;
    int  max_r_new_idx;

    real N;
    real A;
    vector[K] alpha;
        
    sum_w_pot <- 0;

    for (v in 1:N_pot)
      sum_w_pot <- sum_w_pot + d_pot[v,2] * exp(-square(d_pot[v,1])/square(psi));  
    

    w <- exp(-(d2)/square(psi));
  
    for (i in 1:N_cores){    
      // local piece
      for (k in 1:K)
	r_new[i] <- gamma[k] * r[idx_cores[i]];
   
      for (k in 1:K) {out_sum[k] <- 0;}
      sum_w <- 0;
    
      for (k in 1:K){
	for (j in 1:N_cells){ // change N_hood to N_cells
	  if (j != idx_cores[i]){
	    out_sum[k] <- out_sum[k] + w[j,i] * r[j][k];
	  }  
	}
      }
   
      //local vs. non-local
      for (k in 1:K)
	r_new[i,k] <- r_new[i,k] + out_sum[k] * (1-gamma[k]) / sum_w_pot;
    
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
      alpha <- phi .* r_new[i];
      
      A <- sum(alpha);
      N <- sum(y[i]);     

      log_lik[i] <- lgamma(N + 1) + lgamma(A) - lgamma(N + A);
      for (k in 1:K) log_lik[i] <- log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]);
    } 
  }
}