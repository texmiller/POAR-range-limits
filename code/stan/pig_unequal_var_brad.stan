functions {
  // modified Bessel function of the first kind
  /*
    real besselI(real z, real v){
      real sums;
      sums = 0;
      
      for (m in 1:101) {
        if (m + v > 0) {
          sums = sums + pow(0.5 * z, v + 2 * (m - 1)) / ( tgamma(m) * tgamma(m+v) );
        }
        if (v <= -101) {
          reject("v must be > -101 in besselI(z, v)");
        }
      }
      return sums;
    }
  
  // modified Bessel function of the second kind (sometimes referred to as
                                                  // modified Bessel function of the third kind)
  real besselK(real z, real v){
    real v1;
    v1 = v - 1E-12;
    return pi() / 2 * (besselI(z,-v1) - besselI(z,v1)) / sin(v1 * pi());
  }*/
    
    // modified Bessel function of the first kind
    real besselI(real x, real nu) {
    real half_x;
    real log_half_x;
    real initial;
    real summand;
    real biggest;
    real smallest;
    real piece;
    int m;
    int mp1;
    
    half_x = 0.5 * x;
    log_half_x = log(half_x);
    m = 0;
    mp1 = 1;
    initial = 0;
    while ((mp1 + nu) < 0) {
      initial = initial + half_x ^ (2 * m + nu) / 
        (tgamma(mp1) * tgamma(mp1 + nu));
      m = mp1;
      mp1 = m + 1;
    }
    biggest = -lgamma(mp1) -lgamma(mp1 + nu) + (2 * m + nu) * log_half_x;
    m = mp1;
    piece = positive_infinity();
    smallest = -745.13321911;
    summand = 0.0;
    while (piece > smallest) {
      mp1 = m + 1;
      piece = -lgamma(mp1) - lgamma(mp1 + nu) + (2 * m + nu) * log_half_x - biggest;
      summand = summand + exp(piece);
      m = mp1;
    }
    return exp(biggest + log1p(summand)) + initial;
  }
  
  // modified Bessel function of the second kind (sometimes referred to as
                                                  // modified Bessel function of the third kind)
  real besselK(real x, real nu) {
    return 0.5 * pi() * (besselI(x, -nu) - besselI(x, nu)) /
      sin(nu * pi());
  }
  
  // log(pmf) of the Sichel distribution, using the recurrence formula from
  // Chesson and Lee 2005.
  row_vector log_sichel_prob(real omega, real xi, real gamma){
    // Only fitting data from patches 0:50
    row_vector[51] P;
    real alpha;
    real bKd;
    int nn;
    real omega_over_alpha;
    real xi_times_omega_over_alpha;
    real square_xi_times_omega_over_alpha;
    real two_times_xi_times_omega_over_square_alpha;
    
    // calculate constants
    alpha <- sqrt(square(omega) + 2 * omega * xi);
    bKd   <- besselK(omega, gamma);
    omega_over_alpha <- omega / alpha;
    xi_times_omega_over_alpha <- xi * omega_over_alpha;
    square_xi_times_omega_over_alpha <- square(xi_times_omega_over_alpha);
    two_times_xi_times_omega_over_square_alpha <- 2 * xi_times_omega_over_alpha / alpha;
    
    // calculate probabilities
    P[1] <-  1 * pow(omega_over_alpha, gamma)   * besselK(alpha, gamma)   / (tgamma(1) * bKd);
    P[2] <- xi * pow(omega_over_alpha, gamma+1) * besselK(alpha, gamma+1) / (tgamma(2) * bKd);
    for(n in 3:51){
      nn <- n-1;
      P[n] <- two_times_xi_times_omega_over_square_alpha * ( (gamma + nn - 1) / nn ) * P[n-1] + square_xi_times_omega_over_alpha * P[n-2] / ( nn * (nn-1));
    }
    return log(P);
  }
  
    int sichel_rng(row_vector all_log_sichel){
    real r;
    int n;
    row_vector[51] cumsum;
    
    // generate random number
    r <- uniform_rng(0,1);
    
    // calculate cumulative sum of distribution
    cumsum <- cumulative_sum(exp(all_log_sichel));

    // initiate n at 0
    n <- 0;
    
    // find where random number landed
    for(c in 1:51){
       if (r >= cumsum[c]) {
           n <- n+1;
      }
    }
    return n;
  }
}

data {
 int<lower=1> N;
 int<lower=1> L;
 int<lower=1> fit_i[N];
 int<lower=1> colB;
 int<lower=1> colR;
 int<lower=0> y[N];
 matrix[L,colB] MB;
 matrix[L,colR] MR;
 int<lower=1> trt[colR];
}

parameters {
 vector[colB] B_x;
 vector[colB] B_o;
 vector[colR] R_x;
 vector[colR] R_o;
 vector<lower=0>[2] tau_x;
 vector<lower=0>[2] tau_o;
}

transformed parameters {
  vector[L] mu;
  vector[L] sig2;
  vector<lower=0>[L] xi;
  vector<lower=0>[L] omega;
  
  // calculate omega and xi for all combinations of fixed and random effects:
  xi    <- exp(MB*B_x + MR*R_x);
  omega <- exp(MB*B_o + MR*R_o);
  
  // the mean and variance below are only valid when gamma = -0.5
  mu   <- xi;
  sig2 <- xi .* (1 + xi ./ omega);
}

model {
  real fit[N];
  row_vector [51] all_log_sichel[L];
  real gamma;
  
  gamma <- -0.5;
  
  for(i in 1:L){
    all_log_sichel[i] <- log_sichel_prob(omega[i], xi[i], gamma);
  }
  
  tau_o ~ cauchy(0,5);
  tau_x ~ cauchy(0,5);
 
  // R_ values for unequal variance among groups
  for(i in 1:colR){
    R_o[i] ~ normal(0,tau_o[trt[i]]);
    R_x[i] ~ normal(0,tau_x[trt[i]]);
  }

  // R_ values for equal variance among groups
  //R_o ~ normal(0,tau_o);
  //R_x ~ normal(0,tau_x);  
  
  for(i in 1:N) {
    fit[i] <- all_log_sichel[fit_i[i]][y[i] + 1 ];
  }
  
  increment_log_prob(sum(fit));
}

generated quantities{
  real bpvalue;
  row_vector [51] all_log_sichel[L];
  vector[N] log_lik;
  int<lower=0> z[N];
  vector[N] pres_y;
  vector[N] pres_z;
  
  for(i in 1:L){
    all_log_sichel[i] <- log_sichel_prob(omega[i], xi[i], -0.5);
  }
  
  for(i in 1:N) {
    // calculate log-likelihood
    log_lik[i] <- all_log_sichel[fit_i[i]][y[i] + 1 ];
    
    // generate posterior-predictive samples
    z[i]   <- sichel_rng(all_log_sichel[fit_i[i]]);
    
    // calculate pearson residual for real and simulated data
    pres_y[i] <- (y[i] - xi[fit_i[i]]) / sqrt(sig2[fit_i[i]]);
    pres_z[i] <- (z[i] - xi[fit_i[i]]) / sqrt(sig2[fit_i[i]]);
  }
  
  // calculate summed squares for real and simulated data
  bpvalue <- step( sum(pres_y .* pres_y) - sum(pres_z .* pres_z) );
}