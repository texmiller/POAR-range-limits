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
  int<lower=0> n;
  int<lower=0> y[n]; 
}

parameters {
  real<lower=0> lambda; //poisson mean
  real<lower=0> sigma; // IG shape
  real<lower=0> theta[n]; //observation-level deviates
}

model {
  
  // prior
  for(i in 1:n){
    theta[i] ~ ig(1, sigma);
  }

  // sampling
  // need to loop for zero truncation
  for (i in 1:n) {
    y[i] ~ poisson(lambda * theta[i]); //T[1,] for zero truncation
  }
  
}


