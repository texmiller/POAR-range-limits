//this model is for troubleshooting POAR growth
functions {
  
  real ig_lpdf (real x, real mu, real lambda){
    //vector [num_elements (x)] prob;
    real prob;
    real lprob;
    prob = (lambda/(2*pi()*(x^3)))^0.5*exp(-lambda*(x - mu)^2/(2*mu^2*x));
    //for (i in 1:num_elements(x)) {
      //   prob[i] = (lambda/(2*pi()*(x[i]^3)))^0.5*exp(-lambda*(x[i] - mu)^2/(2*mu^2*x[i]));
      //}
    lprob = log( prob ); 
    //lprob = sum (log(prob)); 
    return lprob;
  }
  
}

data {
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations
  int<lower=0> n_g;    // N. of data points for the survival model
  int<lower=0> n_blocks_g;         // N. of blocks
  int<lower=0> site_g[n_g];  // site index
  int<lower=0> block_g[n_g];  // site index
  int<lower=0> source_g[n_g];  // source index
  int<lower=1> y_g[n_g]; // # tillers at time t+1.
  vector[n_g] size_g;  // log size at time t
  vector[n_g] male_g;  // sex (male=1, fem=0)
  vector[n_g] long_g;  // longitude of site
}

parameters {
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
  real bsizelong2_g;
  real blong2sex_g;
  real bsizelong2sex_g;
  //random effects
  real<lower=0> site_tau_g; 
  real site_rfx_g[n_sites];
  real<lower=0> block_tau_g; 
  real block_rfx_g[n_blocks_g];  
  real<lower=0> source_tau_g; 
  real source_rfx_g[n_sources];
  //dispersion parameters
  real d0_g;
  real dsize_g;
  real dsex_g;
  real dlong_g;
  
  real<lower=0> sigma; // IG shape
  real<lower=0> theta[n_g]; //observation-level deviates
}

transformed parameters{
  real predG[n_g];
  real dispG[n_g];
  
  for(igrow in 1:n_g){
    
    // prediction for growth mean
    predG[igrow] = exp( b0_g + 
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
                bsizelong2_g * size_g[igrow] * pow(long_g[igrow],2) +
                blong2sex_g * pow(long_g[igrow],2) * male_g[igrow] +
                bsizelong2sex_g * size_g[igrow] * pow(long_g[igrow],2) * male_g[igrow] +
                //random effects
                site_rfx_g[site_g[igrow]] +
                block_rfx_g[block_g[igrow]] +
                source_rfx_g[source_g[igrow]] );
    // prediction for growth dispersion -- keep this simpler
    dispG[igrow] = exp(d0_g + 
                  dsize_g * size_g[igrow] + 
                  dsex_g * male_g[igrow] + 
                  dlong_g * long_g[igrow]);
  }
  
}

model {
  
  // priors on parameters
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
  bsizelong2_g ~ normal(0, 100);   
  blong2sex_g ~ normal(0, 100);   
  bsizelong2sex_g ~ normal(0, 100);    
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
  d0_g ~ normal(0, 100);
  dsize_g ~ normal(0, 100);
  dsex_g ~ normal(0, 100);
  dlong_g ~ normal(0, 100);

  // prior
  for(i in 1:n_g){
    theta[i] ~ ig(dispG[i], sigma);
  }

  // sampling
  // need to loop for zero truncation
  for (i in 1:n_g) {
    y_g[i] ~ poisson(predG[i] * theta[i]) T[1,]; //T[1,] for zero truncation
  }
  
}


