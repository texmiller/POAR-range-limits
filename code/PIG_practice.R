library(rstan)
library(actuar)

## goal: figure out how to implement the Poisson Inverse Gaussian,
## which has greater skewness and kurtosis than the NB

## first try this user-defined Stan function for the IG, which I found here:
## https://discourse.mc-stan.org/t/writing-a-user-defined-function-for-inverse-gaussian-distribution/1815
## I am a little confused why this person used an integer data set, since the IG is continuous
failure_data <- list (J = 88,
                      y = rep (c (8, 16, 32, 40, 56, 60, 64, 72, 80, 96, 104, 108, 112, 114, 
                                  120, 128, 136, 152, 156, 160, 168, 176, 184, 194, 208, 216, 
                                  224, 232, 240, 246, 256, 264, 272, 280, 288, 304, 308, 328,
                                  340, 352, 358, 360, 384, 392, 400, 424, 438, 448, 464, 480, 
                                  536, 552, 576, 608, 656, 716), times = c (1, 4, 2, 4, 3, 1, 1, 
                                                                            5, 4, 2, 1, 1, 2, 1, 1, 1, 1, 3, 1, 1, 5, 1, 3, 1, 2, 1, 4, 1,
                                                                            1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                                            1, 1, 1, 1, 1, 1, 1)) 
)
#pdf of inverse Gaussian distribution
#f(x; mu, shape) = (shape/(2*pi*x^3))^0.5*exp(-shape*(x - mu)^2/(2*mu^2*x))

model_IG <- "functions {

  real IG_log (vector x, real mu, real lambda){
vector [num_elements (x)] prob;
real lprob;
for (i in 1:num_elements(x)) {
  prob[i] <- (lambda/(2*pi()*(x[i]^3)))^0.5*exp(-lambda*(x[i] - mu)^2/(2*mu^2*x[i]));
}
lprob <- sum (log(prob)); 
return lprob;

}
}

data {
int<lower=0> J; // number of observations
vector[J] y; // failure in hours

}
parameters {
real<lower=0> mu; // mean
real<lower=0> lambda; // shape

}

model {
//define priors of mu and lambda


//Likelihood of the data
y ~ IG(mu, lambda);
}
"

mod_IG <- stan_model(model_code = model_IG)

fit_IG <- sampling (mod_IG, data = failure_data, chains = 4, iter = 10^4,
                    warmup = 5000, thin = 20)
fit_IG

## ok this works! Can I make it a PIG?
## simulate data from PIG distribution 
pig_dat <- list(y = rpoisinvgauss(100, 
                                  mean=10, 
                                  shape = 0.2),
                N = 100)

## try fitting in Stan as Poisson-IG mixture, using the IG function from above
model_PIG <- "functions {

  real IG_log (vector x, real mu, real shape){
vector [num_elements (x)] prob;
real lprob;
for (i in 1:num_elements(x)) {
  prob[i] <- (shape/(2*pi()*(x[i]^3)))^0.5*exp(-shape*(x[i] - mu)^2/(2*mu^2*x[i]));
}
lprob <- sum (log(prob)); 
return lprob;

}
}

data {
int<lower=0> N; // number of observations
int<lower=0> y[N]; // failure in hours
}
parameters {
real<lower=0> lambda; //poisson mean
real<lower=0> sigma; // IG shape
vector[N] theta; //observation-level deviates
}

model {
theta ~ IG(1, sigma);
for(i in 1:N){
y[i] ~ poisson(lambda*theta[i]);
}

}
"
mod_PIG <- stan_model(model_code = model_PIG)
fit_PIG <- sampling (mod_PIG, data = pig_dat, chains = 3, iter = 10^4,
                    warmup = 5000, thin = 20)
## this does not work
## I get this error indicating that the poisson mean is negative,
## which suggests the theta's are negative, and I don't see how this can be
## since the IG has support (0,Inf)