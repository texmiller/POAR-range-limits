
data { 

  // int<lower=0> n_sites;  // N. of locations
  int<lower=0> n_blocks; // N. of locations
  int<lower=0> n_s;      // N. of data points for the survival model  
  
  int<lower=0> site_s[n_s];  // site placeholder
  int<lower=0> block_s[n_s]; // site placeholder
  
  // Data for the survival model
  vector[n_s] size_s; // log size at time t
  vector[n_s] male_s; // male or not?
  vector[n_s] long_s; // longitude
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  
}

parameters {
  
  real block_a_u_s;
  real<lower=0> block_a_tau_s;
  vector[n_blocks] block_a_s_zzz;
  
  real b_s;   // Survival reg. slope
  real b_l;   // Survival reg. slope
  real b_sex;   // Survival intercept for MALES
  real b_sex_l; // Male by longitude interaction
  
}

transformed parameters {
  
  // placeholders of year intercept random year effects
  vector[n_blocks] block_a_s;
  real mS[n_s];
  
  block_a_s = block_a_u_s + block_a_tau_s * block_a_s_zzz;
  
  // prediction
  for(nsurv in 1:n_s){
    mS[nsurv] = block_a_s[block_s[nsurv]] + 
                b_s   * size_s[nsurv] +
                b_l   * long_s[nsurv] +
                b_sex * male_s[nsurv] +
                b_sex_l * male_s[nsurv] * long_s[nsurv];
  }
  
}

model {
  
  // Hyper-priors
  block_a_u_s   ~ normal(0, 100);
  block_a_s_zzz ~ normal(0, 10);
  block_a_tau_s ~ inv_gamma(0.001, 0.001);
  
  // block ranef
  for(j in 1:n_blocks){
    block_a_s_zzz[j] ~ normal(block_a_u_s, block_a_tau_s);
  }
  
  // priors on parameters
  b_s     ~ normal(0, 100);   // Survival slope
  b_l     ~ normal(0, 100);   // Longitude slope
  b_sex   ~ normal(0, 100);   // Sex slope
  b_sex_l ~ normal(0, 100);   // Sex slope
 
  // sampling
  y_s ~ bernoulli_logit(mS);
  
}

generated quantities{
    
  vector[n_s] log_lik;
  
  for(i in 1:n_s){
    log_lik[i] = bernoulli_logit_lpmf(y_s[i] | mS[i]);
  }
    
}
