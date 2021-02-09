data {
  int<lower=0> N;         // number of pumps
  int y[N];             // number of failiures
  real t[N];     // operation time of pump
}
parameters {
  real<lower=0> lambda[N];       // expected count

  real<lower=0> alpha0;             // population treatment effect
  real<lower=0> beta0;              // population treatment effect
}
model {
  target+=exponential_lpdf(alpha0 | 1.0);
  target+=gamma_lpdf(beta0 | 0.1, 1.0);
  for (i in 1:N){
    target+=gamma_lpdf(lambda[i] | alpha0, beta0);
    target+=poisson_lpmf(y[i] | lambda[i]*t[i]);
  }
}

