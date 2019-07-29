data{
  
  int<lower=0> n_g;  
  
  int<lower=0> y_germ[n_g];  // site placeholder
  int<lower=0> tot_seeds[n_g]; // site placeholder
  
  // Data for the survival model
  real SR[n_g]; // log size at time t
  
}

parameters{
  
  real<lower=0.00001,upper=0.99999> v0;
  real<lower=0.00001> a_g;
  
}

transformed parameters{
  
  real<lower=0,upper=1> prob_g[n_g];
  
  for(ii in 1:n_g){
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

