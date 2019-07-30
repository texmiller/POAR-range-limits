
data { 

  int<lower=0> n_sites;         // N. of sites
  
  // Data for survival model
  int<lower=0> n_blocks_s;         // N. of locations
  int<lower=0> n_s;    // N. of data points for the survival model  
  
  int<lower=0> site_s[n_s];  // site index
  int<lower=0> block_s[n_s]; // block indexes index
  
  vector[n_s] size_s;  // log size at time t
  vector[n_s] male_s;  // male or not?
  vector[n_s] long_s;  // male or not?
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.

  // Data for growth model
  int<lower=0> n_blocks_g;         // N. of locations
  int<lower=0> n_g;    // N. of data points for the growth model  
  
  int<lower=0> site_g[n_g];  // site index
  int<lower=0> block_g[n_g]; // block indexes index
  
  vector[n_g] size_g;  // log n. of tillers at time t
  vector[n_g] male_g;  // male or not?
  vector[n_g] long_g;  // longitude (Centered)
  int<lower=0> y_g[n_g]; // n. of tillers at time t+1.
  
  // Data for flowering probability model
  int<lower=0> n_f;         // N. of data points for flowering mod.  
  int<lower=0> n_blocks_f;  // N. of blocks
  
  int<lower=0> site_f[n_f];  // site index
  int<lower=0> block_f[n_f]; // block indexes
  
  vector[n_f] size_f;  // log size at time t+1
  vector[n_f] male_f;  // male or not?
  vector[n_f] long_f;  // longitude (centered)
  int<lower=0,upper=1> y_f[n_f]; // Flowering at time t+1.
  
  // Data for panicules model
  int<lower=0> n_blocks_p;  // N. of blocks
  int<lower=0> n_p;         // N. of data points for the panicule model  
  
  int<lower=0> site_p[n_p];  // site index
  int<lower=0> block_p[n_p]; // block indexes index
  
  vector[n_p] size_p;  // log n. of tillers at time t+1
  vector[n_p] male_p;  // male or not?
  vector[n_p] long_p;  // longitude (Centered)
  int<lower=0> y_p[n_p]; // n. of panicules at time t+1.
  
  // Data for seed viability model
  int<lower=0> n_v;   // data points
  int<lower=0> y_v[n_v];  // site placeholder
  int<lower=0> tot_seeds[n_v]; // site placeholder
  real SR[n_v]; // Sex ratio
  
}

parameters {
  
  // Survival parameters
  real site_a_u_s;
  real<lower=0> site_a_tau_s; 
  real site_a_s[n_sites];
  
  real b_s_s;   // Survival size slope
  real b_l_s;   // Survival longitude slope
  real b_sex_s;   // Survival intercept for MALES
  real b_sex_l_s; // Male by longitude interaction
  
  // Growth parameters
  real site_a_u_g;
  real<lower=0> site_a_tau_g; 
  real site_a_g[n_sites];
  
  real b_s_g;   // Growth size slope
  real b_l_g;   // Growth longitude slope
  real b_sex_g;   // Growth intercept for MALES
  real b_sex_l_g; // Male by longitude interaction
  real<lower=0> phi_g; // Growth dispersion parameters 
  
  // Flowering probability parameters
  real site_a_u_f;
  real<lower=0> site_a_tau_f; 
  real site_a_f[n_sites];
  
  real b_s_f;   // Growth size slope
  real b_l_f;   // Growth longitude slope
  real b_sex_f;   // Growth intercept for MALES
  real b_sex_l_f; // Male by longitude interaction
  
   // Panicules parameters
  real site_a_u_p;
  real<lower=0> site_a_tau_p; 
  real site_a_p[n_sites];
  
  real b_s_p;   // Growth size slope
  real b_l_p;   // Growth longitude slope
  real b_sex_p;   // Growth intercept for MALES
  real b_sex_l_p; // Male by longitude interaction
  real<lower=0> phi_p; // Growth dispersion parameters 
  
  // Seed viability parameters
  real<lower=0,upper=1> v0;
  real<lower=0> a_v;
  
}

transformed parameters{
  
  real pred_s[n_s];
  real pred_g[n_g];
  real pred_f[n_f];
  real pred_p[n_p];
  real<lower=0,upper=1> pred_v[n_v];
  
  // prediction for survival
  for(is in 1:n_s){
    pred_s[is] = site_a_s[site_s[is]] + 
                 b_s_s     * size_s[is] +
                 b_l_s     * long_s[is] +
                 b_sex_s   * male_s[is] +
                 b_sex_l_s * male_s[is] * long_s[is];
  }
  
  // prediction for growth
  for(ig in 1:n_g){
    pred_g[ig] = site_a_g[site_g[ig]] + 
                 b_s_g     * size_g[ig] +
                 b_l_g     * long_g[ig] +
                 b_sex_g   * male_g[ig] +
                 b_sex_l_g * male_g[ig] * long_g[ig];
  }
  
  // prediction for flowering probability
  for(iif in 1:n_f){
    pred_f[iif] = site_a_f[site_f[iif]] + 
                  b_s_f     * size_f[iif] +
                  b_l_f     * long_f[iif] +
                  b_sex_f   * male_f[iif] +
                  b_sex_l_f * male_f[iif] * long_f[iif];
  }
 
  // prediction for growth
  for(ip in 1:n_p){
    pred_p[ip] = site_a_p[site_p[ip]] + 
                 b_s_p     * size_p[ip] +
                 b_l_p     * long_p[ip] +
                 b_sex_p   * male_p[ip] +
                 b_sex_l_p * male_p[ip] * long_p[ip];
  } 
  
  // Prediction for seed viability
  for(ii in 1:n_v){
    pred_v[ii] = v0 * (1 - pow(SR[ii],a_v) ) + 0.00001;
  }
  
}

model {
  
  // Hyper-priors (0 because slopes should be ~0)
  site_a_u_s    ~ normal(0, 100);
  site_a_tau_s  ~ inv_gamma(0.001, 0.001);
  site_a_u_g    ~ normal(0, 100);
  site_a_tau_g  ~ inv_gamma(0.001, 0.001);
  site_a_u_f    ~ normal(0, 100);
  site_a_tau_f  ~ inv_gamma(0.001, 0.001);
  site_a_u_p    ~ normal(0, 100);
  site_a_tau_p  ~ inv_gamma(0.001, 0.001);

  // hyperparameters years
  for (i in 1:n_sites){
    site_a_s[i] ~ normal(site_a_u_s, site_a_tau_s);
    site_a_g[i] ~ normal(site_a_u_g, site_a_tau_g);
    site_a_f[i] ~ normal(site_a_u_f, site_a_tau_f);
    site_a_p[i] ~ normal(site_a_u_p, site_a_tau_p);
  }

  // priors on parameters
  b_s_s     ~ normal(0, 100);   // Survival slope
  b_l_s     ~ normal(0, 100);   // Longitude slope
  b_sex_s   ~ normal(0, 100);   // Sex slope
  b_sex_l_s ~ normal(0, 100);   // Sex slope
  
  b_s_g     ~ normal(0, 100);   // growth size slope
  b_l_g     ~ normal(0, 100);   // Longitude slope
  b_sex_g   ~ normal(0, 100);   // Sex slope (for growth)
  b_sex_l_g ~ normal(0, 100);   // Sex by long. slope (growth)

  b_s_f     ~ normal(0, 100);   // Flowering size slope
  b_l_f     ~ normal(0, 100);   // Longitude slope (flow.)
  b_sex_f   ~ normal(0, 100);   // Sex slope (for flow.)
  b_sex_l_f ~ normal(0, 100);   // Sex by long. slope (for flow.)

  b_s_p     ~ normal(0, 100);   // panicules size slope
  b_l_p     ~ normal(0, 100);   // Longitude slope
  b_sex_p   ~ normal(0, 100);   // Sex slope (for panicules)
  b_sex_l_p ~ normal(0, 100);   // Sex by long. slope (panicules)

  v0  ~ beta(10,1);  // intercept viability model
  a_v ~ gamma(1,1);  // "decay" in viability with SR

  // Sampling 
  y_s ~ bernoulli_logit(pred_s);
  y_g ~ neg_binomial_2_log(pred_g, phi_g);
  y_f ~ bernoulli_logit(pred_f);
  y_p ~ neg_binomial_2_log(pred_p, phi_p);
  y_v ~ binomial(tot_seeds, pred_v);
  
}
