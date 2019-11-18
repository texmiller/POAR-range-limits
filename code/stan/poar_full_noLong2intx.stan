
data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations
  
  // Data for survival sub-model (s)
  int<lower=0> n_s;    // N. of data points for the survival model
  int<lower=0> n_blocks_s;         // N. of blocks
  int<lower=0> site_s[n_s];  // site index
  int<lower=0> block_s[n_s];  // site index
  int<lower=0> source_s[n_s];  // source index
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  vector[n_s] size_s;  // log size at time t
  vector[n_s] male_s;  // sex (male=1, fem=0)
  vector[n_s] long_s;  // longitude of site
  
  // Data for growth sub-model (g)
  int<lower=0> n_g;    // N. of data points for the survival model
  int<lower=0> n_blocks_g;         // N. of blocks
  int<lower=0> site_g[n_g];  // site index
  int<lower=0> block_g[n_g];  // site index
  int<lower=0> source_g[n_g];  // source index
  int<lower=1> y_g[n_g]; // # tillers at time t+1.
  vector[n_g] size_g;  // log size at time t
  vector[n_g] male_g;  // sex (male=1, fem=0)
  vector[n_g] long_g;  // longitude of site
    
  // Data for flowering sub-model (f)
  int<lower=0> n_f;    // N. of data points for the survival model
  int<lower=0> n_blocks_f;         // N. of blocks
  int<lower=0> site_f[n_f];  // site index
  int<lower=0> block_f[n_f];  // site index
  int<lower=0> source_f[n_f];  // source index
  int<lower=0,upper=1> y_f[n_f]; // Flowering at time t+1.
  vector[n_f] size_f;  // log size at time t
  vector[n_f] male_f;  // sex (male=1, fem=0)
  vector[n_f] long_f;  // longitude of site
    
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
  
  // Data for seed viability sub-model (v)
  int<lower=0> n_v;   // data points
  int<lower=0> y_v[n_v];  // number of viable seeds
  int<lower=0> tot_seeds_v[n_v]; // number of trials
  real SR_v[n_v]; // Sex ratio (proportion female?)
  // Data for seed germination sub-model (m)
  int<lower=0> n_m;   // data points
  int<lower=0> y_m[n_m];  // number of germinating seeds
  int<lower=0> tot_seeds_m[n_m]; // number of trials
  real SR_m[n_m]; // Sex ratio (proportion female?)
  
  // data for seed count
  int<lower=0> n_d;   // data points
  int<lower=0> y_d[n_d];  // total number of seeds per panicle 
}

parameters {
  //Survival
  //fixed effects
  real b0_s;   
  real bsize_s;   
  real bsex_s;   
  real blong_s;  
  real bsizesex_s;  
  real bsizelong_s;
  real blongsex_s;
  real bsizelongsex_s;
  real blong2_s;  
  //real bsizelong2_s;
  //real blong2sex_s;
  //real bsizelong2sex_s;
  //random effects
  real<lower=0> site_tau_s; 
  real site_rfx_s[n_sites];
  real<lower=0> block_tau_s; 
  real block_rfx_s[n_blocks_s];  
  real<lower=0> source_tau_s; 
  real source_rfx_s[n_sources];
  
  //Growth
  //fixed effects
  real b0_g;   
  real bsize_g;  
  real bsex_g;   
  real blong_g;   
  real bsizesex_g;   
  real bsizelong_g;
  real blongsex_g;
  real bsizelongsex_g;
  real blong2_g;  
  //real bsizelong2_g;
  //real blong2sex_g;
  //real bsizelong2sex_g;
  //random effects
  real<lower=0> site_tau_g; 
  real site_rfx_g[n_sites];
  real<lower=0> block_tau_g; 
  real block_rfx_g[n_blocks_g];  
  real<lower=0> source_tau_g; 
  real source_rfx_g[n_sources];
  real<lower=0> phi_g; // Growth dispersion parameter 
    
  //Flowering
  //fixed effects
  real b0_f;    // Survival size intercept
  real bsize_f;   // Survival size slope
  real bsex_f;   // Survival sex effect
  real blong_f;   // Survival longitude slope
  real bsizesex_f;   // Interactions
  real bsizelong_f;
  real blongsex_f;
  real bsizelongsex_f;
  real blong2_f;  
  //real bsizelong2_f;
  //real blong2sex_f;
  //real bsizelong2sex_f;
  //random effects
  real<lower=0> site_tau_f; 
  real site_rfx_f[n_sites];
  real<lower=0> block_tau_f; 
  real block_rfx_f[n_blocks_f];  
  real<lower=0> source_tau_f; 
  real source_rfx_f[n_sources];
      
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
  //real bsizelong2_p;
  //real blong2sex_p;
  //real bsizelong2sex_p;
  //random effects
  real<lower=0> site_tau_p; 
  real site_rfx_p[n_sites];
  real<lower=0> block_tau_p; 
  real block_rfx_p[n_blocks_p];  
  real<lower=0> source_tau_p; 
  real source_rfx_p[n_sources];
  real<lower=0> phi_p; // Panicle dispersion parameter
  
  // Seed viability parameters
  real<lower=0,upper=1> v0;
  real<lower=0> a_v;
  real<lower=0.1> phi_v;             
  // Germination rate
  real<lower=0,upper=1> m;
  real<lower=0.1> phi_m;  
  // Seed count
  real<lower=0> lambda_d;
  }

transformed parameters{
  
  real predS[n_s];
  real predG[n_g];
  real predF[n_f];
  real predP[n_p];
  real<lower=0,upper=1> predV[n_v];
  real<lower=0,upper=1> predM[n_m];
  
  // beta-binom reparameterization for viab and germ
  vector[n_v] alpha_v; 
  vector[n_v] beta_v;   
  vector[n_m] alpha_m; 
  vector[n_m] beta_m;  
  
  // prediction for survival
  for(isurv in 1:n_s){
    predS[isurv] = b0_s + 
                //main effects
                bsize_s * size_s[isurv] + bsex_s * male_s[isurv] + blong_s * long_s[isurv] + 
                //2-way interactions
                bsizesex_s * size_s[isurv] * male_s[isurv] +
                bsizelong_s * size_s[isurv] * long_s[isurv] +
                blongsex_s * long_s[isurv] * male_s[isurv] +
                //3-way interaction
                bsizelongsex_s * size_s[isurv] * long_s[isurv] * male_s[isurv] +
                //quadratic longitude effects
                blong2_s * pow(long_s[isurv],2) + 
                //bsizelong2_s * size_s[isurv] * pow(long_s[isurv],2) +
                //blong2sex_s * pow(long_s[isurv],2) * male_s[isurv] +
                //bsizelong2sex_s * size_s[isurv] * pow(long_s[isurv],2) * male_s[isurv] +
                //random effects
                site_rfx_s[site_s[isurv]] +
                block_rfx_s[block_s[isurv]] +
                source_rfx_s[source_s[isurv]];
  }
  
  // prediction for growth
  for(igrow in 1:n_g){
    predG[igrow] = b0_g + 
                //main effects
                bsize_g * size_g[igrow] + bsex_g * male_g[igrow] + blong_g * long_g[igrow] + 
                //2-way interactions
                bsizesex_g * size_g[igrow] * male_g[igrow] +
                bsizelong_g * size_g[igrow] * long_g[igrow] +
                blongsex_g * long_g[igrow] * male_g[igrow] +
                //3-way interaction
                bsizelongsex_g * size_g[igrow] * long_g[igrow] * male_g[igrow] +
                //quadratic longitude effects
                blong2_g * pow(long_g[igrow],2) + 
                //bsizelong2_g * size_g[igrow] * pow(long_g[igrow],2) +
                //blong2sex_g * pow(long_g[igrow],2) * male_g[igrow] +
                //bsizelong2sex_g * size_g[igrow] * pow(long_g[igrow],2) * male_g[igrow] +
                //random effects
                site_rfx_g[site_g[igrow]] +
                block_rfx_g[block_g[igrow]] +
                source_rfx_g[source_g[igrow]];
  }

  // prediction for flowering
  for(iflow in 1:n_f){
    predF[iflow] = b0_f + 
                //main effects
                bsize_f * size_f[iflow] + bsex_f * male_f[iflow] + blong_f * long_f[iflow] + 
                //2-way interactions
                bsizesex_f * size_f[iflow] * male_f[iflow] +
                bsizelong_f * size_f[iflow] * long_f[iflow] +
                blongsex_f * long_f[iflow] * male_f[iflow] +
                //3-way interaction
                bsizelongsex_f * size_f[iflow] * long_f[iflow] * male_f[iflow] +
                //quadratic longitude effects
                blong2_f * pow(long_f[iflow],2) + 
                //bsizelong2_f * size_f[iflow] * pow(long_f[iflow],2) +
                //blong2sex_f * pow(long_f[iflow],2) * male_f[iflow] +
                //bsizelong2sex_f * size_f[iflow] * pow(long_f[iflow],2) * male_f[iflow] +
                //random effects
                site_rfx_f[site_f[iflow]] +
                block_rfx_f[block_f[iflow]] +
                source_rfx_f[source_f[iflow]];
  }

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
                //bsizelong2_p * size_p[ipan] * pow(long_p[ipan],2) +
                //blong2sex_p * pow(long_p[ipan],2) * male_p[ipan] +
                //bsizelong2sex_p * size_p[ipan] * pow(long_p[ipan],2) * male_p[ipan] +
                //random effects
                site_rfx_p[site_p[ipan]] +
                block_rfx_p[block_p[ipan]] +
                source_rfx_p[source_p[ipan]];
  }  
  
  // Prediction for seed viability
  for(iviab in 1:n_v){
    predV[iviab] = v0 * (1 - pow(SR_v[iviab],a_v) ) + 0.00001;
    alpha_v[iviab] = predV[iviab] * phi_v;
    beta_v[iviab] = (1 - predV[iviab]) * phi_v; 
    }
  // Prediction for germination
  for(igerm in 1:n_m){
    predM[igerm] = m * v0 * (1 - pow(SR_m[igerm],a_v) ) + 0.00001;
    alpha_m[igerm] = predM[igerm] * phi_m;
    beta_m[igerm] = (1 - predM[igerm]) * phi_m; 
    }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // priors on parameters
  //Survival
  b0_s ~ normal(0, 100);   
  bsize_s ~ normal(0, 100);   
  bsex_s ~ normal(0, 100);   
  blong_s ~ normal(0, 100);   
  bsizesex_s ~ normal(0, 100);   
  bsizelong_s ~ normal(0, 100);   
  blongsex_s ~ normal(0, 100);   
  bsizelongsex_s ~ normal(0, 100);   
  blong2_s ~ normal(0, 100);   
  //bsizelong2_s ~ normal(0, 100);   
  //blong2sex_s ~ normal(0, 100);   
  //bsizelong2sex_s ~ normal(0, 100);    
  site_tau_s ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sites){
    site_rfx_s[i] ~ normal(0, site_tau_s);
  }
  block_tau_s ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_blocks_s){
    block_rfx_s[i] ~ normal(0, block_tau_s);
  }
  source_tau_s ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sources){
    source_rfx_s[i] ~ normal(0, source_tau_s);
  }
  //Growth
  b0_g ~ normal(0, 100);   
  bsize_g ~ normal(0, 100);   
  bsex_g ~ normal(0, 100);   
  blong_g ~ normal(0, 100);   
  bsizesex_g ~ normal(0, 100);   
  bsizelong_g ~ normal(0, 100);   
  blongsex_g ~ normal(0, 100);   
  bsizelongsex_g ~ normal(0, 100);   
  blong2_g ~ normal(0, 100);   
  //bsizelong2_g ~ normal(0, 100);   
  //blong2sex_g ~ normal(0, 100);   
  //bsizelong2sex_g ~ normal(0, 100);    
  site_tau_g ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sites){
    site_rfx_g[i] ~ normal(0, site_tau_g);
  }
  block_tau_g ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_blocks_g){
    block_rfx_g[i] ~ normal(0, block_tau_g);
  }
  source_tau_g ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sources){
    source_rfx_g[i] ~ normal(0, source_tau_g);
  }
  //Flowering
  b0_f ~ normal(0, 100);   
  bsize_f ~ normal(0, 100);   
  bsex_f ~ normal(0, 100);   
  blong_f ~ normal(0, 100);   
  bsizesex_f ~ normal(0, 100);   
  bsizelong_f ~ normal(0, 100);   
  blongsex_f ~ normal(0, 100);   
  bsizelongsex_f ~ normal(0, 100);   
  blong2_f ~ normal(0, 100);   
  //bsizelong2_f ~ normal(0, 100);   
  //blong2sex_f ~ normal(0, 100);   
  //bsizelong2sex_f ~ normal(0, 100);    
  site_tau_f ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sites){
    site_rfx_f[i] ~ normal(0, site_tau_f);
  }
  block_tau_f ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_blocks_f){
    block_rfx_f[i] ~ normal(0, block_tau_f);
  }
  source_tau_f ~ inv_gamma(0.001, 0.001);
  for (i in 1:n_sources){
    source_rfx_f[i] ~ normal(0, source_tau_f);
  }
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
  //bsizelong2_p ~ normal(0, 100);   
  //blong2sex_p ~ normal(0, 100);   
  //bsizelong2sex_p ~ normal(0, 100);  
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
  
  v0  ~ beta(10,1);  // intercept viability model
  a_v ~ gamma(1,1);  // "decay" in viability with SR
  phi_v ~ pareto(0.1,1.5);
  m  ~ beta(10,1);  // intercept viability model
  phi_m ~ pareto(0.1,1.5);
  lambda_d ~ inv_gamma(0.001, 0.001);
  
  // sampling
  for (i in 1:n_s) {
  y_s[i] ~ bernoulli_logit(predS[i]);
  }
  for (i in 1:n_g) {
  y_g[i] ~ neg_binomial_2_log(predG[i], phi_g);
  target += - log1m(neg_binomial_2_log_lpmf(0 | predG[i], phi_g)); // manually zero-truncating
  }
  for (i in 1:n_f) {
  y_f[i] ~ bernoulli_logit(predF[i]);
  }
  for (i in 1:n_p) {
  y_p[i] ~ neg_binomial_2_log(predP[i], phi_p);
  target += - log1m(neg_binomial_2_log_lpmf(0 | predP[i], phi_p)); // manually zero-truncating
  }
  for (i in 1:n_v) {
  y_v[i] ~ beta_binomial(tot_seeds_v[i], alpha_v[i], beta_v[i]);
  }
  for (i in 1:n_m) {
  y_m[i] ~ beta_binomial(tot_seeds_m[i], alpha_m[i], beta_m[i]);
  }
  for (i in 1:n_d){
  y_d[i] ~ poisson(lambda_d);  
  }

  
}


