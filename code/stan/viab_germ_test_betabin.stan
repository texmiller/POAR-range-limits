

// The input data is a vector 'y' of length 'N'.
data {
  // Data for seed viability sub-model (v)
  int<lower=0> n_v;   // data points
  int<lower=0> y_v[n_v];  // number of viable seeds
  int<lower=0> tot_seeds_v[n_v]; // number of trials
  real SR_v[n_v]; // Sex ratio (proportion female?)
  // Data for seed germination sub-model (g)
  int<lower=0> n_g;   // data points
  int<lower=0> y_g[n_g];  // number of germinating seeds
  int<lower=0> tot_seeds_g[n_g]; // number of trials
  real SR_g[n_g]; // Sex ratio (proportion female?)

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // Seed viability parameters
  real<lower=0,upper=1> v0;
  real<lower=0> a_v;
  real<lower=0.1> phi_v;             
  // Germination rate
  real<lower=0,upper=1> g;
  real<lower=0.1> phi_g;             
}

transformed parameters{
  real<lower=0,upper=1> predV[n_v];
  real<lower=0,upper=1> predG[n_g];
  
  // beta-binom reparameterization
  vector[n_v] alpha_v; 
  vector[n_v] beta_v;   
  vector[n_g] alpha_g; 
  vector[n_g] beta_g;  
  
  // Prediction for seed viability
  for(iviab in 1:n_v){
    predV[iviab] = v0 * (1 - pow(SR_v[iviab],a_v) ) + 0.00001;
    alpha_v[iviab] = predV[iviab] * phi_v;
    beta_v[iviab] = (1 - predV[iviab]) * phi_v;  
  }
  // Prediction for germination
  for(igerm in 1:n_g){
    predG[igerm] = g * v0 * (1 - pow(SR_g[igerm],a_v) ) + 0.00001;
    alpha_g[igerm] = predG[igerm] * phi_g;
    beta_g[igerm] = (1 - predG[igerm]) * phi_g;  
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  v0  ~ beta(10,1);  // intercept viability model
  a_v ~ gamma(1,1);  // "decay" in viability with SR
  phi_v ~ pareto(0.1,1.5);
  g  ~ beta(10,1);  // intercept viability model
  phi_g ~ pareto(0.1,1.5);
  
  y_v ~ beta_binomial(tot_seeds_v, alpha_v, beta_v);
  y_g ~ beta_binomial(tot_seeds_g, alpha_g, beta_g);
}

