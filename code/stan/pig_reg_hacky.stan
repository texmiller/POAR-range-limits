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
  int<lower=0> N; // number of observations
  int<lower=0> y[N]; // failure in hours
  vector[N] x; // failure in hours
}
parameters {
  real b0; 
  real b1; 
  real<lower=0> sigma; // IG shape
  real<lower=0> theta[N]; //observation-level deviates
}

transformed parameters{
  vector<lower=0>[N] yhat;
  
  yhat = exp(b0 + b1 * x );
  
}

model{
  real mu;
  mu = 1;
  
  // prior
  for(i in 1:N){
    theta[i] ~ ig(mu, sigma);
  }
  
  for(i in 1:N){
    y[i] ~ poisson( yhat[i] * theta[i] );
  }
  
}
