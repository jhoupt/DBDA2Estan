data {
  int<lower=0> n_total;
  int<lower=1>n_x1_lvl;
  int<lower=1>n_x2_lvl;
  real y[n_total];
  int x1[n_total];
  int x2[n_total];
  real<lower=0> a_gamma_sh_ra[2];
}
transformed data { 
  real y_mean;
  real y_sd;
  y_mean = mean(y);
  y_sd = sd(y);
}
parameters { 
  real<lower=0> y_sigma;
  real a0;
  vector[n_x1_lvl] a1;
  vector[n_x2_lvl] a2;
  vector[n_x2_lvl] a1a2[n_x1_lvl];
  real<lower=0> a1_sd;
  real<lower=0> a2_sd;
  real<lower=0> a1a2_sd;
}
transformed parameters { 
  vector[n_total] mu;
  for (n in 1:n_total) { 
    mu[n] = a0 + a1[x1[n]] + a2[x2[n]] + a1a2[x1[n],x2[n]];
  }
}
model {
  y ~ normal(mu, y_sigma);
  //y_sigma ~ uniform(y_sd / 100, y_sd * 10);

  a0 ~ normal(y_mean, y_sd * 5);

  a1 ~ normal(0.0, a1_sd);
  a2 ~ normal(0.0, a2_sd);
  for (n in 1:n_x1_lvl) { 
    a1a2[n] ~ normal(0.0, a1a2_sd);
  }

  // or try a folded t (Cauchy)
  a1_sd ~ gamma(a_gamma_sh_ra[1], a_gamma_sh_ra[2]);
  a2_sd ~ gamma(a_gamma_sh_ra[1], a_gamma_sh_ra[2]);
  a1a2_sd ~ gamma(a_gamma_sh_ra[1], a_gamma_sh_ra[2]);
}

generated quantities {
  matrix[n_x1_lvl,n_x2_lvl] m;
  real b0;
  vector[n_x1_lvl] b1;
  vector[n_x2_lvl] b2;
  real b1b2[n_x1_lvl,n_x2_lvl];

  // Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
  for (j1 in 1:n_x1_lvl) { 
    for (j2 in 1:n_x2_lvl) {
      m[j1,j2] = a0 + a1[j1] + a2[j2] + a1a2[j1,j2]; // cell means 
    } 
  }
  b0 = mean(m);
  for (j1 in 1:n_x1_lvl) { b1[j1] = mean(row(m,j1)) - b0; }
  for (j2 in 1:n_x2_lvl) { b2[j2] = mean(col(m,j2)) - b0; }
  for (j1 in 1:n_x1_lvl) { 
    for (j2 in 1:n_x2_lvl) {
      b1b2[j1,j2] = m[j1,j2] - (b0 + b1[j1] + b2[j2]);
    } 
  }    
}

