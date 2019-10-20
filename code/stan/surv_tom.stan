//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_blocks;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations
  int<lower=0> n_s;    // N. of data points for the survival model
  int<lower=0> site_s[n_s];  // site index
  int<lower=0> block_s[n_s];  // site index
  int<lower=0> source_s[n_s];  // source index
  
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  vector[n_s] size_s;  // log size at time t
  vector[n_s] male_s;  // sex (male=1, fem=0)
  vector[n_s] long_s;  // longitude of site
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //fixed effects
  real b_0;    // Survival size intercept
  real b_size;   // Survival size slope
  real b_sex;   // Survival sex effect
  real b_long;   // Survival longitude slope
  real b_size_sex;   // Interactions
  real b_size_long;
  real b_long_sex;
  real b_size_long_sex;

  //random effects
  real<lower=0> site_tau; 
  real site_rfx[n_sites];
  real<lower=0> block_tau; 
  real block_rfx[n_blocks];  
  real<lower=0> source_tau; 
  real source_rfx[n_sources];
  }

transformed parameters{
  
  real mS[n_s];
  // prediction for survival
  for(isurv in 1:n_s){
    mS[isurv] = b_0 + 
                #main effects
                b_size * size_s[isurv] + b_sex * male_s[isurv] + b_long * long_s[isurv] + 
                #2-way interactions
                b_size_sex * size_s[isurv] * male_s[isurv] +
                b_size_long * size_s[isurv] * long_s[isurv] +
                b_long_sex * long_s[isurv] * male_s[isurv] +
                #3-way interaction
                b_size_long_sex * size_s[isurv] * long_s[isurv] * male_s[isurv] +
                #random effects
                site_rfx[site_s[isurv]] +
                block_rfx[block_s[isurv]] +
                source_rfx[source_s[isurv]];
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // priors on parameters
  b_0 ~ normal(0, 100);   
  b_size ~ normal(0, 100);   
  b_sex ~ normal(0, 100);   
  b_long ~ normal(0, 100);   
  b_size_sex ~ normal(0, 100);   
  b_size_long ~ normal(0, 100);   
  b_long_sex ~ normal(0, 100);   
  b_size_long_sex ~ normal(0, 100);  
  
  site_tau ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sites){
    site_rfx[i] ~ normal(0, site_tau);
  }
  block_tau ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_blocks){
    block_rfx[i] ~ normal(0, block_tau);
  }
  source_tau ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sources){
    source_rfx[i] ~ normal(0, source_tau);
  }
  
  // sampling
  y_s ~ bernoulli_logit(mS);
}

generated quantities {
 real y_s_new[n_s];

 for (i in 1:n_s) {
 y_s_new[i] = bernoulli_logit_rng(mS[i]);
 }

}
