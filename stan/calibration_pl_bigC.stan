// Author:  Andria Dawson
// Date:    26 June 2014
// Settlement era pollen estimation model based on STEPPS1
// Uses veg proportions and pollen counts to estimate process parameters
// With a Gaussian dispersal model using a single dispersal distance parameter psi
 

data {
  int<lower=0> K;                // number of taxa
  int<lower=0> N_cores;          // number of core sites
  int<lower=0> N_cells;          // number of spatial cells
  int<lower=0> N_pot;           // number of potential contributing spatial cells

  int y[N_cores,K];              // pollen count data
  
  int idx_cores[N_cores];        // core cell indices
  
  vector[K] r[N_cells];          // pls proportions
  //  matrix[N_cells,N_cores] d2;    // distance matrix squared
  matrix[N_cells,N_cores] d;    // distance matrix squared
  matrix[N_pot, 2] d_pot;        // distances and counts of potential neighborhood
}

transformed data {
}

parameters {
  vector<lower=0.01, upper=300>[K] phi;  // dirichlet precision pars
  // real<lower=0.1, upper=2> psi;          // dispersal par
  real<lower=0, upper=1> gamma;          // local proportion par
  real<lower=0, upper=500> a;
  real<lower=2, upper=100> b;
}

transformed parameters {

}

model {

  // declarations
  matrix[N_cells,N_cores] w; 
  vector[N_pot] w_pot;       
  vector[K] r_new[N_cores];
  vector[K] out_sum;    
  
  real sum_w;
  real C;
  real max_r_new;
  int  max_r_new_idx;

  //print("psi = ", psi);

  // priors 
  phi      ~ uniform(0.01,300);
  //psi      ~ uniform(0.1,2);
  gamma    ~ uniform(0,1);
  a ~ uniform(0, 500);   
  b ~ uniform(2, 100) ; 

  C <- 0;
  for (v in 1:N_pot)
    C <- C + d_pot[v,2] * (b-2) * (b-1) / ( 2 * pi() * a * a ) * pow( 1 + d_pot[v,1] / a, -b) ;//exp(-square(d_pot[v,1])/square(psi));
  //print("C = ", C);  
  
  //  w <- exp(-(d2)/square(psi));
  for (i in 1:N_cells)
    for (j in 1:N_cores)
      w[i,j] <- (b-2) * (b-1) / ( 2 * pi() * a * a ) * pow( 1 + d[i,j] / a, -b) ;
      
  for (i in 1:N_cores){   
   /*
    if (sum(r[idx_cores[i]]) == 0 ){
      print("Cell ", idx_cores[i], " has NO vegetation! Neighborhood props will not sum to 1.");
    }
    */
    
    // local piece
    r_new[i] <- gamma*r[idx_cores[i]];
   
    for (k in 1:K) {out_sum[k] <- 0;}
    sum_w <- 0;
    for (j in 1:N_cells){ // change N_hood to N_cells
      if (j != idx_cores[i]){
        out_sum <- out_sum + w[j,i]*r[j];
	sum_w   <- sum_w + w[j,i];
       }  
     }

    //print("C = ", C);
    //print("sum_w = ", sum_w);
     
    //sum_w   <- sum(out_sum);
  
    //local vs. non-local
    //r_new[i] <- r_new[i] + out_sum*(1-gamma)/sum_w;
    r_new[i] <- r_new[i] + out_sum * (1-gamma) / C;
    
    //print(r_new[i]);
    //print(sum(r_new[i]));
    

    //print("sum r_new ", sum(r_new[i]));
 
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
