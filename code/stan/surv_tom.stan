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
  int<lower=0> n_s;    // N. of data points for the survival model
  vector[n_s] size_s;  // log size at time t
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b_0;    // Survival size intercept
  real b_s;   // Survival size slope

}

transformed parameters{
  
  real mS[n_s];
  // prediction for survival
  for(nsurv in 1:n_s){
    mS[nsurv] = b_0 + b_s * size_s[nsurv];
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y_s ~ bernoulli_logit(mS);
}

