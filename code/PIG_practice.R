library(rstan)
library(actuar)
library(dplyr)
library(statmod)
library(gamlss)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )


# user defined exponential -----------------------------------------------

stan_d_l <- list( N = 1000, 
                  y = rexp(1000, 2.577) )

model_exp <- 
  
  "functions {
  
    real expp_lpdf (vector x, real lambda){
    
      vector [num_elements (x)] prob;
      real lprob;
      for (i in 1:num_elements(x)) {
        prob[i] = log(lambda) - x[i] * lambda;
      }
      lprob = sum( prob );
      return lprob;
  
    }
  
  }
  
  data {
    int<lower=0> N; // number of observations
    vector[N] y; // failure in hours
  }
  parameters {
    real<lower=0> lambda; // rate
  }
  
  model {
    //define priors of mu and lambda
    
    //Likelihood of the data
    y ~ expp(lambda);
  }
  "

mod_expp <- stan_model(model_code = model_exp)
fit_expp <- sampling( mod_expp, data = stan_d_l, 
                      chains = 4, iter = 10^4,
                      warmup = 5000, thin = 20 )


# Inverse Gaussian -------------------------------------------------

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


model_ig <- 
  "functions {
  
    real ig_lpdf (vector x, real mu, real lambda){
    
      vector [num_elements (x)] prob;
      real lprob;
      for (i in 1:num_elements(x)) {
        prob[i] = (lambda/(2*pi()*(x[i]^3)))^0.5*exp(-lambda*(x[i] - mu)^2/(2*mu^2*x[i]));
      }
      lprob = sum( log(prob) );
      return lprob;
  
    }
  
  }
  
  data {
    int<lower=0> N; // number of observations
    vector[N] y; // failure in hours
  }
  
  parameters {
    real<lower=0> mu; // mean
    real<lower=0> lambda; // shape
  }
  
  model {
    //Likelihood of the data
    y ~ ig(mu, lambda);
  }
  "

n_rep  <- 1000
ig_sim <- list( N = n_rep,
                y = rinvgauss(n_rep, 100, 50) )
mod_ig <- stan_model(model_code = model_ig)
fit_ig <- sampling (mod_ig, data = ig_sim, chains = 4, iter = 10^4,
                    warmup = 5000, thin = 20)
fit_IG


# Poisson-Inverse Gaussian (pig) -----------------------------------------------

## ok this works! Can I make it a PIG?
## simulate data from PIG distribution 

n_rep   <- 1000
pig_dat <- list(y = rpoisinvgauss(1000, 
                                  mean=10, 
                                  shape = 0.2),
                N = 1000)

## try fitting in Stan as Poisson-IG mixture, using the IG function from above
model_pig <- 
"functions {

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
}
parameters {
  real<lower=0> lambda; //poisson mean
  real<lower=0> sigma; // IG shape
  real<lower=0> theta; // 
  // vector<lower=0>[N] theta; //observation-level deviates
}

model{
  // prior
  theta ~ ig(1, sigma);
  
  for(i in 1:N){
    y[i] ~ poisson( lambda * theta );
  }

}
"

mod_PIG <- stan_model(model_code = model_pig)
fit_PIG <- sampling( mod_PIG, 
                     data   = pig_dat, 
                     chains = 3, 
                     iter   = 10^4,
                     warmup = 5000, 
                     thin   = 20 )


# PIG parameterization from Ben Goodrich ------------------------


pig_goodrich <- 
'functions {
  // ignoring the 2pi constant
  real ig_lpdf(real x, real mu, real shape){
    return 0.5 * log(shape) - 1.5 * log(x) - shape * square( (x - mu) / mu) / x;
  }
}

data {
  int<lower=1> N;
  int<lower=0> x[N];
}

parameters {
  real<lower=0> mu;
  real<lower=0> shape;
  real<lower=0> lambda;
}

model {
  x ~ poisson(lambda);
  lambda ~ ig(mu, shape); // edited. Thanks Bob!
}'

# generate inverse Gamma issue
n_rep   <- 100
pig_dat <- list(x = rpoisinvgauss(n_rep, 
                                  mean=10, 
                                  shape = 0.2),
                N = n_rep)

# model pig from Goodrich
mod_pig_gr <- stan_model(model_code = pig_goodrich)
fit_PIG <- sampling( mod_pig_gr, 
                     data   = pig_dat, 
                     chains = 3, 
                     iter   = 10^4,
                     warmup = 5000, 
                     thin   = 20 )



# Poisson-Inverse Gaussian (pig) regression -----------------------------------------------

# this fits, but I don't think it makes sense
n_rep   <- 33
x       <- rep( seq(1,30,by=1), 33)
pig_dat <- list(y = rpoisinvgauss( length(x), 
                                   mean  = exp(x*0.2) , 
                                   shape = 20),
                x = x,
                N = length(x) )

plot(pig_dat$x, pig_dat$y)


model_pig_reg <- 
"functions {

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
  int<lower=0>  N; // number of observations
  int<lower=0>  y[N]; // n. tillers t1
  real<lower=0> x[N]; // n. tillers t0
}
parameters {
  //real<lower=0> lambda; //poisson mean
  real<lower=0> sigma; // IG shape
  real<lower=0> theta[N]; // 
  real beta; // 
  // vector<lower=0>[N] theta; //observation-level deviates
}

transformed parameters{
  real<lower=0> lambda[N]; //poisson mean
  
  for(i in 1:N)
    lambda[i] = exp(beta * x[i]);
  
}

model{
  // prior
  for(i in 1:N){
    theta[i] ~ ig(lambda[i], sigma);
  }
  
  for(i in 1:N){
    y[i] ~ poisson( lambda[i] * theta[i] );
  }

}
"

mod_pig_r <- stan_model(model_code = model_pig_reg)
fit_PIG   <- sampling( mod_pig_r, 
                       data   = pig_dat, 
                       chains = 3, 
                       iter   = 4000,
                       warmup = 1000, 
                       thin   = 3 )

plot( summary(fit_PIG)$summary[paste0('lambda[',1:990,']'),][,'mean'],
      pig_dat$x )



# PIG model of the MEAN (using Poission Generalized Ingerse Gamma (PGIG)) -----

# simulate data
N       <- 200
mu      <- 3
shape   <- 2  # changing shape to sigma to match gamlss's parameterization
sigma   <- sqrt(1/shape)
x       <- replace( x, x == 0, 1 )
y       <- rPIG( N, 
                 mu    = 10, 
                 sigma = 2 )
pig_dat <- list( y     = y,
                 y_max = max(y) + 1,
                 N     = N )

# compile model
mod_pig   <- stan_model( file = 'code/stan/pgig_univ.stan' )
# mod_pig   <- readRDS( file = 'code/stan/pgig_univ.rds' )

# fit model
fit_PIG   <- sampling( object = mod_pig, 
                       data   = pig_dat, 
                       chains = 1, 
                       iter   = 4000,
                       warmup = 1000, 
                       thin   = 3,
                       verbose = T)

fit_PIG



# PIG estimation of GROUP effects (via PGIG) --------------------------------


# simulate data
N       <- 100
x       <- runif(N, 1, 5) %>% round(0) # set up groups
y       <- rPIG( N, mu = x, sigma = 1 ) # response variable
pig_dat <- list( y     = y,
                 x     = x,
                 x_i   = x,
                 x_u   = unique(x) %>% sort,
                 x_n   = n_distinct(x),
                 y_max = max(y) + 1,
                 N     = N )

# compile model
mod_pig_r  <- stan_model( file = 'code/stan/pgig_groups.stan' )
# mod_pig_r  <- readRDS( file = 'code/stan/pgig_groups.rds' )

# fit model
fit_PIG    <- sampling(  
                        mod_pig_r, 
                        data    = pig_dat, 
                        chains  = 1, 
                        iter    = 4000,
                        warmup  = 1000, 
                        thin    = 2,
                        verbose = T 
                      )

fit_PIG


# PIG regression via PGIG --------------------------------------------

# simulate data
N       <- 130
x       <- runif(N, 1, 5) %>% round(0)
y       <- rPIG( N, 
                 mu = exp(x * 0.3), # exponential function 
                 sigma = 1 )
pig_dat <- list( y     = y,
                 x     = x,
                 x_i   = x,
                 x_u   = unique(x) %>% sort,
                 x_n   = n_distinct(x),
                 y_max = max(y) + 1,
                 N     = N )

# does data look like an actual regression?
plot( jitter(pig_dat$x), pig_dat$y )

# compile model
mod_pig_r  <- stan_model( file = 'code/stan/pgig_reg.stan' )
mod_pig_r  <- readRDS( file = 'code/stan/pgig_reg.rds' )

# fit model
fit_PIG    <- sampling( mod_pig_r, 
                        data    = pig_dat, 
                        chains  = 1, 
                        iter    = 4000,
                        warmup  = 1000, 
                        thin    = 2,
                        verbose = T 
                        )

fit_PIG
