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
}

parameters {
  vector<lower=0.01, upper=300>[K] phi;  // dirichlet precision pars
  real<lower=0, upper=1> gamma;          // local proportion par

  vector<lower=0.00001, upper=500>[K] a;
  vector<lower=2, upper=100>[K] b;
}

transformed parameters {
}

model {

  // declarations
  matrix[N_cells,N_cores] w[K];      
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  real sum_w;
  vector[K] sum_w_pot;
  real max_r_new;
  int  max_r_new_idx;

  // priors 
  phi     ~ uniform(0.01,300);
  gamma   ~ uniform(0,1); 

  for (k in 1:K){
      a[k] ~ uniform(0.00001, 500);   
      b[k] ~ uniform(2, 100) ;
  }  
  
  for (k in 1:K){
    sum_w_pot[k] <- 0;
  
    for (v in 1:N_pot)
      sum_w_pot <- sum_w_pot + d_pot[v,2] * (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] ) * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
  }

  for (k in 1:K){
    for (i in 1:N_cells)
      for (j in 1:N_cores)
	w[k][i,j] <- (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] ) * pow( 1 + d[i,j] / a[k], -b[k]) ;
  }
  
  for (i in 1:N_cores){    
    // local piece
    r_new[i] <- gamma*r[idx_cores[i]];
   
    for (k in 1:K) {out_sum[k] <- 0;}
    sum_w <- 0;
    
    for (k in 1:K){
      for (j in 1:N_cells){ // change N_hood to N_cells
	if (j != idx_cores[i]){
	  out_sum[k] <- out_sum[k] + w[k][j,i]*r[j][k];
	}  
      }
    }

    //local vs. non-local
    for (k in 1:K)
      r_new[i,k] <- r_new[i,k] + out_sum[k]*(1-gamma)/sum_w_pot[k];
    
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
// generated quantities {
//   vector[N_cores] log_lik;

//   {
//   // declarations
//   matrix[N_cells,N_cores] w[K];      
//   vector[K] r_new[N_cores];
//   vector[K] out_sum;    
  
//   real sum_w;
//   vector[K] sum_w_pot;
//   real max_r_new;
//   int  max_r_new_idx;

//   real N;
//   real A;
//   vector[K] alpha;  
  
//   for (k in 1:K){
//     sum_w_pot[k] <- 0;
  
//     for (v in 1:N_pot)
//       sum_w_pot <- sum_w_pot + d_pot[v,2] * (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] ) * pow( 1 + d_pot[v,1] / a[k], -b[k]) ; 
//   }

//   for (k in 1:K){
//     for (i in 1:N_cells)
//       for (j in 1:N_cores)
// 	w[k][i,j] <- (b[k]-2) * (b[k]-1) / ( 2 * pi() * a[k] * a[k] ) * pow( 1 + d[i,j] / a[k], -b[k]) ;
//   }
  
//   for (i in 1:N_cores){
//     /*   for (j in 1:N_hood){
// 	 if (idx_hood[i,j] > 0){
// 	 w[i,j] <- exp(-square(d[i,idx_hood[i,j]]/psi));
// 	 } else {
// 	 w[i,j] <- 0;
// 	 }
// 	 }
//     */ 
   
//     /*
//       if (sum(r[idx_cores[i]]) == 0 ){
//       print("Cell ", idx_cores[i], " has NO vegetation! Neighborhood props will not sum to 1.");
//       }
//     */
    
//     // local piece
//     r_new[i] <- gamma*r[idx_cores[i]];
   
//     for (k in 1:K) {out_sum[k] <- 0;}
//     sum_w <- 0;
    
//     for (k in 1:K){
//       for (j in 1:N_cells){ // change N_hood to N_cells
// 	if (j != idx_cores[i]){
// 	  out_sum[k] <- out_sum[k] + w[k][j,i]*r[j][k];
// 	}  
//       }
//     }

//     //sum_w   <- sum(out_sum);
     
//     //print("out_sum", sum(out_sum));
//     //print("sum_w", sum_w);
    
//     //local vs. non-local
//     for (k in 1:K)
//       r_new[i,k] <- r_new[i,k] + out_sum[k]*(1-gamma)/sum_w_pot[k];
    
//     //print(r_new[i]);
//     //print(sum(r_new[i]));
    
//     /*
//       if (sum(r_new[i]) < 1-1e-6 ){
//       print("i ", i);
//       print("core props ", r[idx_cores[i]]);
//       print("sum r_new ", sum(r_new[i]));
//       }
//     */
    
//     // hacky!
//     // find taxon with highest proportional value
//     max_r_new <- 0;
//     for (k in 1:K){
//       if (r_new[i,k] > max_r_new){
//         max_r_new     <- r_new[i,k];
//         max_r_new_idx <- k;
//       }
//     }
    
//     for (k in 1:K){
//       if (r_new[i,k] == 0){
//         r_new[i,k] <- 0.0001;
//         r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
//         print("warning: zero proportion; core: ", i, "; taxon: ", k, " -> adjusting");
//       }
//     }
   
//     alpha <- phi .* r_new[i];
      
//     A <- sum(alpha);
//     N <- sum(y[i]);     

//     log_lik[i] <- lgamma(N + 1) + lgamma(A) - lgamma(N + A);
//     for (k in 1:K) log_lik[i] <- log_lik[i] - lgamma(y[i,k] + 1) + lgamma(y[i,k] + alpha[k]) - lgamma(alpha[k]);
//   }
    
// } 
// }