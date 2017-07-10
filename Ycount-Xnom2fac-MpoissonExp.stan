data {
  int<lower=1> n_cell;
  int<lower=0> y[n_cell];
  int x1[n_cell];
  int x2[n_cell];
  int<lower=1>n_x1_lvl;
  int<lower=1>n_x2_lvl;
  real y_log_mean;
  real y_log_sd;
  real a_gamma_sh_ra[2];
}
parameters {
  real a0;
  vector[n_x1_lvl] a1;
  vector[n_x2_lvl] a2;
  vector[n_x2_lvl] a1a2[n_x1_lvl];
  real<lower=0> a1_sd;
  real<lower=0> a2_sd;
  real<lower=0> a1a2_sd;
}
transformed parameters { 
  vector<lower=0>[n_cell] lambda;
  for (n in 1:n_cell) { 
    lambda[n] = a0 + a1[x1[n]] + a2[x2[n]] + a1a2[x1[n],x2[n]];
  }
  lambda = exp(lambda);
}
model {
  y ~ poisson(lambda);
  a0 ~ normal(y_log_mean, y_log_sd*2);
  
  a1 ~ normal(0.0, a1_sd);
  a1_sd ~ gamma(a_gamma_sh_ra[1],a_gamma_sh_ra[2]);
  
  a2 ~ normal(0.0, a2_sd);
  a2_sd ~ gamma(a_gamma_sh_ra[1],a_gamma_sh_ra[2]);
  
  for (n in 1:n_x1_lvl) { 
    a1a2[n] ~ normal(0.0, a1a2_sd);
  }
  a1a2_sd ~ gamma(a_gamma_sh_ra[1],a_gamma_sh_ra[2]);
}
generated quantities { 
  matrix[n_x1_lvl,n_x2_lvl] m;
  real b0;
  vector[n_x1_lvl] b1;
  vector[n_x2_lvl] b2;
  real b1b2[n_x1_lvl,n_x2_lvl];
  matrix[n_x1_lvl,n_x2_lvl] expm;
  matrix[n_x1_lvl,n_x2_lvl] pp_x1x2_p;
  vector<lower=0>[n_x1_lvl] pp_x1_p;
  vector<lower=0>[n_x2_lvl] pp_x2_p;

  // Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
  for (j1 in 1:n_x1_lvl) { 
    for (j2 in 1:n_x2_lvl) {
      m[j1,j2] = a0 + a1[j1] + a2[j2] + a1a2[j1,j2]; // cell means 
    } 
  }

  b0 = mean(m);
  for (j1 in 1:n_x1_lvl) { b1[j1] = mean(row(m,j1)) - b0; }
  for (j2 in 1:n_x2_lvl) { b2[j2] = mean(col(m,j2)) - b0; }
  for (j1 in 1:n_x1_lvl) { for (j2 in 1:n_x2_lvl) {
    b1b2[j1,j2] = m[j1,j2] - (b0 + b1[j1] + b2[j2]);
  } }    

  // Compute predicted proportions:
  expm = exp(m);
  pp_x1x2_p = expm/sum(expm);

  for (j1 in 1:n_x1_lvl) { pp_x1_p[j1] = sum(row(pp_x1x2_p,j1)); }
  for (j2 in 1:n_x2_lvl) { pp_x2_p[j2] = sum(col(pp_x1x2_p,j2)); }
}
