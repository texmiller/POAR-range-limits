data{
  
  int<lower=0> n_v;  
  
  int<lower=0> y_germ[n_v];  // site placeholder
  int<lower=0> tot_seeds[n_v]; // site placeholder
  
  real SR[n_v]; // Sex ratio
  
}

parameters{
  
  real<lower=0,upper=1> v0;
  real<lower=0> a_g;
  
}

transformed parameters{
  
  real<lower=0,upper=1> prob_g[n_v];
  
  for(ii in 1:n_v){
    prob_g[ii] = v0 * (1 - pow(SR[ii],a_g) ) + 0.00001;
  }
  
}

model{
  
  // priors
  v0  ~ beta(10,1);
  a_g ~ gamma(1,1);
  
  // sampling
  y_germ ~ binomial(tot_seeds, prob_g);
  
}

