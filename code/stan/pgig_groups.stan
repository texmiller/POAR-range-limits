functions {
  // modified Bessel function of the first kind
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
  }
  
  // log(pmf) of the Sichel distribution, using the recurrence formula from
  // Chesson and Lee 2005.
  vector log_sichel_prob(real omega, real xi, real gamma, int y_max){
    // Only fitting data from patches 0:50
    vector[y_max] P;
    real alpha;
    real bKd;
    int nn;
    real omega_over_alpha;
    real xi_times_omega_over_alpha;
    real square_xi_times_omega_over_alpha;
    real two_times_xi_times_omega_over_square_alpha;
    
    // calculate constants
    alpha = sqrt(square(omega) + 2 * omega * xi);
    bKd   = besselK(omega, gamma);
    omega_over_alpha = omega / alpha;
    xi_times_omega_over_alpha = xi * omega_over_alpha;
    square_xi_times_omega_over_alpha = square(xi_times_omega_over_alpha);
    two_times_xi_times_omega_over_square_alpha = 2 * xi_times_omega_over_alpha / alpha;
    
    // calculate probabilities
    P[1] =  1 * pow(omega_over_alpha, gamma)   * besselK(alpha, gamma)   / (tgamma(1) * bKd);
    P[2] = xi * pow(omega_over_alpha, gamma+1) * besselK(alpha, gamma+1) / (tgamma(2) * bKd);
    for(n in 3:y_max){
      nn = n-1;
      P[n] = two_times_xi_times_omega_over_square_alpha * ( (gamma + nn - 1) / nn ) * P[n-1] + square_xi_times_omega_over_alpha * P[n-2] / ( nn * (nn-1));
    }
    return log(P);
  }
}

data {
  int<lower=1> N;
  int<lower=0> y[N];
  int<lower=0> y_max;
  int<lower=0> x_i[N];
  int<lower=0> x_n;
}

parameters {
  real<lower=0, upper=2> omega;
  vector<lower=0, upper=10>[x_n] xi;
  //real<lower=0, upper=10> gamma;
}

transformed parameters {
  vector[x_n] mu;
  vector[x_n] sig2;
  
  // the mean and variance below are only valid when gamma = -0.5
  mu   = xi;
  sig2 = xi .* (1+ xi ./ omega);
}

model {
  real fit[N];
  vector[y_max] logsichel;
  
  for(i in 1:N) {
    logsichel = log_sichel_prob(omega, xi[x_i[i]], -0.5, y_max);
    // y[i] + 1 because 0 corresponds to P[1] in log_sichel_prob
    fit[i] = logsichel[ y[i] + 1 ];
  }
  
  target += sum(fit);
  
}

