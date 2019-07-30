
data { 

  int<lower=0> n_sites;         // N. of locations
  int<lower=0> n_blocks;         // N. of locations
  int<lower=0> n_s;    // N. of data points for the survival model  
  
  int<lower=0> site_s[n_s];  // site index
  int<lower=0> block_s[n_s]; // block indexes index
  
  // Data for the survival model
  vector[n_s] size_s;  // log size at time t
  vector[n_s] male_s;  // male or not?
  vector[n_s] long_s;  // male or not?
  int<lower=0> y_s[n_s]; // Survival at time t+1.
  
}

parameters {
  
  real site_a_u_s;
  real<lower=0> site_a_tau_s; 
  real site_a_s[n_sites];
  
  real b_s;   // Survival reg. slope
  real b_l;   // Survival reg. slope
  real b_sex;   // Survival intercept for MALES
  real b_sex_l; // Male by longitude interaction
  real<lower=0> phi_g; // Growth dispersion parameter
  
}

transformed parameters{
  
  real mS[n_s];
  
  // prediction for survival
  for(nsurv in 1:n_s){
    mS[nsurv] = site_a_s[ site_s[nsurv] ] + 
                b_s   * size_s[nsurv] +
                b_l   * long_s[nsurv] +
                b_sex * male_s[nsurv] +
                b_sex_l * male_s[nsurv] * long_s[nsurv];
  }
  
}

model {
  
  // Hyper-priors (0 because slopes should be ~0)
  site_a_u_s    ~ normal(0, 100);
  site_a_tau_s  ~ inv_gamma(0.001, 0.001);
  
  // priors
  phi_g ~ inv_gamma(0.001, 0.001);  
  
  // hyperparameters years
  for (i in 1:n_sites){
    site_a_s[i] ~ normal(site_a_u_s, site_a_tau_s);
  }

  // priors on parameters
  b_s     ~ normal(0, 100);   // Survival slope
  b_l     ~ normal(0, 100);   // Longitude slope
  b_sex   ~ normal(0, 100);   // Sex slope
  b_sex_l ~ normal(0, 100);   // Sex slope
 
  // Sampling 
  y_s ~ neg_binomial_2_log(mS, phi_g);
  
}

generated quantities{
    
  vector[n_s] log_lik;
  
  for(i in 1:n_s){
    log_lik[i] = neg_binomial_2_log_lpmf(y_s[i] | mS[i], phi_g);
  }
    
}

