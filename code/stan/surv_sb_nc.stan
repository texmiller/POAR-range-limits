
data { 

  int<lower=0> n_sites;  // N. of locations
  int<lower=0> n_blocks; // N. of locations
  int<lower=0> n_s;      // N. of data points for the survival model  
  
  int<lower=0> site_s[n_s];  // site placeholder
  int<lower=0> block_s[n_s]; // site placeholder
  
  // Data for the survival model
  vector[n_s] size_s; // log size at time t
  vector[n_s] male_s; // male or not?
  vector[n_s] long_s; // male or not?
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  
}

parameters {
  
  real site_a_u_s;
  real<lower=0> site_a_tau_s;
  vector[n_sites] site_a_s_zzz;

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
  vector[n_sites] site_a_s;
  vector[n_blocks] block_a_s;
  
  site_a_s  = site_a_u_s  + site_a_tau_s  * site_a_s_zzz;
  block_a_s = block_a_u_s + block_a_tau_s * block_a_s_zzz;
  
}

model {
  
  real mS[n_s];
  int i_s_site;
  int i_s_block;
  
  // Hyper-priors
  site_a_u_s    ~ normal(0, 100);
  site_a_tau_s  ~ inv_gamma(0.001, 0.001);
  block_a_u_s    ~ normal(0, 100);
  block_a_tau_s  ~ inv_gamma(0.001, 0.001);

  // site ranef
  for (i in 1:n_sites){
    site_a_s_zzz[i] ~ normal(site_a_u_s, site_a_tau_s);
  }
  
  for (i in 1:n_blocks){
    block_a_s_zzz[i] ~ normal(block_a_u_s, block_a_tau_s);
  }

  // site ranef
  
  // priors on parameters
  b_s     ~ normal(0, 100);   // Survival slope
  b_l     ~ normal(0, 100);   // Longitude slope
  b_sex   ~ normal(0, 100);   // Sex slope
  b_sex_l ~ normal(0, 100);   // Sex slope
 
  // Sampling -------------------------------------------
  
  // survival
  for(nsurv in 1:n_s){
    i_s_site  = site_s[nsurv];
    i_s_block = block_s[nsurv];
    mS[nsurv] = site_a_s[i_s_site] +
                block_a_s[i_s_block] + 
                b_s   * size_s[nsurv] +
                b_l   * long_s[nsurv] +
                b_sex * male_s[nsurv] +
                b_sex_l * male_s[nsurv] * long_s[nsurv];
  }
  y_s ~ bernoulli_logit(mS);
  
}
