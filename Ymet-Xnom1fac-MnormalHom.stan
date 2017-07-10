data {
  int<lower=1> n_total;
  int x[n_total];
  int y[n_total]; 
  int<lower=1> n_x_lvl;
  real y_mean;
  real<lower=0> y_sd;
}
parameters {
  real a0;
  vector[n_x_lvl] a;
  real<lower=0> a_sigma;
  real<lower=0> y_sigma;
}
model {
  y ~ normal(a0 + a[x], y_sigma^2);
  y_sigma ~ uniform(y_sd/100, y_sd*10);
  a0 ~ normal(y_mean, (y_sd*5)^2);
  #
  a ~ normal(0.0, a_sigma^2)
  a_sigma ~ gamma(a_gamma_sh_ra[1], a_gamma_sh_ra[2]);
}
generated quantities { 
  real b0;
  vector[n_x_lvl] b;
  vector[n_x_lvl] m;

  # Convert a0,a[] to sum-to-zero b0,b[] :
  m = a0 + a;
  b0 = mean(m);
  b = m - b0 ;
}
