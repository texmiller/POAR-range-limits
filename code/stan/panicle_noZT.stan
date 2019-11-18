
data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations

  // Data for fertility sub-model (p)
  int<lower=0> n_p;    // N. of data points for the survival model
  int<lower=0> n_blocks_p;         // N. of blocks
  int<lower=0> site_p[n_p];  // site index
  int<lower=0> block_p[n_p];  // site index
  int<lower=0> source_p[n_p];  // source index
  int<lower=1> y_p[n_p]; // # panicles at time t.
  vector[n_p] size_p;  // log size at time t
  vector[n_p] male_p;  // sex (male=1, fem=0)
  vector[n_p] long_p;  // longitude of site

}

parameters {
  //Panicles
  //fixed effects
  real b0_p;    // Survival size intercept
  real bsize_p;   // Survival size slope
  real bsex_p;   // Survival sex effect
  real blong_p;   // Survival longitude slope
  real bsizesex_p;   // Interactions
  real bsizelong_p;
  real blongsex_p;
  real bsizelongsex_p;
  real blong2_p;  
  real bsizelong2_p;
  real blong2sex_p;
  real bsizelong2sex_p;
  //random effects
  real<lower=0> site_tau_p; 
  real site_rfx_p[n_sites];
  real<lower=0> block_tau_p; 
  real block_rfx_p[n_blocks_p];  
  real<lower=0> source_tau_p; 
  real source_rfx_p[n_sources];
  real<lower=0> phi_p; // Panicle dispersion parameter
  }

transformed parameters{
  
  real predP[n_p];
  
  // prediction for panicles
  for(ipan in 1:n_p){
    predP[ipan] = b0_p + 
                //main effects
                bsize_p * size_p[ipan] + bsex_p * male_p[ipan] + blong_p * long_p[ipan] + 
                //2-way interactions
                bsizesex_p * size_p[ipan] * male_p[ipan] +
                bsizelong_p * size_p[ipan] * long_p[ipan] +
                blongsex_p * long_p[ipan] * male_p[ipan] +
                //3-way interaction
                bsizelongsex_p * size_p[ipan] * long_p[ipan] * male_p[ipan] +
                //quadratic longitude effects
                blong2_p * pow(long_p[ipan],2) + 
                bsizelong2_p * size_p[ipan] * pow(long_p[ipan],2) +
                blong2sex_p * pow(long_p[ipan],2) * male_p[ipan] +
                bsizelong2sex_p * size_p[ipan] * pow(long_p[ipan],2) * male_p[ipan] +
                //random effects
                site_rfx_p[site_p[ipan]] +
                block_rfx_p[block_p[ipan]] +
                source_rfx_p[source_p[ipan]];
  }  
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
    //Panicles
  b0_p ~ normal(0, 100);   
  bsize_p ~ normal(0, 100);   
  bsex_p ~ normal(0, 100);   
  blong_p ~ normal(0, 100);   
  bsizesex_p ~ normal(0, 100);   
  bsizelong_p ~ normal(0, 100);   
  blongsex_p ~ normal(0, 100);   
  bsizelongsex_p ~ normal(0, 100);   
  blong2_p ~ normal(0, 100);   
  bsizelong2_p ~ normal(0, 100);   
  blong2sex_p ~ normal(0, 100);   
  bsizelong2sex_p ~ normal(0, 100);  
  site_tau_p ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sites){
    site_rfx_p[i] ~ normal(0, site_tau_p);
  }
  block_tau_p ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_blocks_p){
    block_rfx_p[i] ~ normal(0, block_tau_p);
  }
  source_tau_p ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sources){
    source_rfx_p[i] ~ normal(0, source_tau_p);
  }
  
  
  // sampling
  //for (i in 1:n_s) {
  //y_s[i] ~ bernoulli_logit(predS[i]);
  //}
  //for (i in 1:n_g) {
  //y_g[i] ~ neg_binomial_2_log(predG[i], phi_g);
  //target += - log1m(neg_binomial_2_log_lpmf(0 | predG[i], phi_g)); // manually zero-truncating
  //}
  //for (i in 1:n_f) {
  //y_f[i] ~ bernoulli_logit(predF[i]);
  //}
  //for (i in 1:n_p) {
  //y_p[i] ~ neg_binomial_2_log(predP[i], phi_p);
  //target += - log1m(neg_binomial_2_log_lpmf(0 | predP[i], phi_p)); // manually zero-truncating
  //}
  //for (i in 1:n_v) {
  //y_v[i] ~ beta_binomial(tot_seeds_v[i], alpha_v[i], beta_v[i]);
  //}
  //for (i in 1:n_m) {
  //y_m[i] ~ beta_binomial(tot_seeds_m[i], alpha_m[i], beta_m[i]);
  //}
  //for (i in 1:n_d){
  //y_d[i] ~ poisson(lambda_d);  
  //}

  y_p ~ neg_binomial_2_log(predP, phi_p);

  
}


