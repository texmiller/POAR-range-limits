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
 int x_i[N]; // index of of predictors (from 1 leaf, to max num. leaves)
 int x_n;    // unique(x_i)
 real x_u[x_n];
}

parameters {
  real b_x;
  real<lower=0,upper=3> omega;
}

transformed parameters {
  vector[x_n] mu;
  vector[x_n] sig2;
  vector[x_n] xi;
  
  // calculate omega and xi for all combinations of fixed and random effects:
  for( i in 1:x_n){
    mu[i]   = exp( b_x * x_u[i] );
  }
  
  xi   = mu;
  // the mean and variance below are only valid when gamma = -0.5
  sig2 = mu .* (1 + mu);
  
}

model {
  real fit[N];
  vector[y_max] logsichel;
  
  // Priors: tightest possible
  b_x ~ normal(0.3, 10);
  omega ~ gamma(1, 1);
  
  for(i in 1:N) {
    logsichel = log_sichel_prob(omega, xi[x_i[i]], -0.5, y_max);
    // y[i] + 1 because 0 corresponds to P[1] in log_sichel_prob
    fit[i] = logsichel[ y[i] + 1 ];
  }
  
  target += sum(fit);
  
}
