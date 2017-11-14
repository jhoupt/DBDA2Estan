data {
  int<lower=1> n_x;
  int<lower=1> n_total;
  vector[n_total] y;
  matrix[n_total,n_x] x;
}
transformed data {
  // Standardize the data:
  matrix[n_total,n_x] zx;
  vector[n_total] zy;
  vector[n_x] x_m;
  vector[n_x] x_sd;
  real y_m;
  real y_sd;

  y_m = mean(y);
  y_sd = sd(y);
  zy = (y - y_m) / y_sd; 

  for (j in 1:n_x) { 
    x_m[j] = mean(col(x,j));
    x_sd[j] = sd(col(x,j));
    zx[,j] = (col(x,j) - x_m[j]) / x_sd[j]; 
  }
}
parameters {
  real z_beta_0;
  vector[n_x] z_beta;
  real<lower=0> z_sigma;
  real<lower=0> sigma_beta;
  real<lower=0> nu;
}

model {

  // Specify the model for the standardized data
  zy ~ student_t(nu, z_beta_0 + zx * z_beta , z_sigma);

  // Priors vague on standardized scale:
  z_beta_0 ~ normal(0, 4);
  z_beta ~ student_t(1, 0, sigma_beta);
  z_sigma ~ uniform(1.0E-5, 1.0E+1);
  nu ~ exponential(1/30.0);
  sigma_beta ~ gamma(2.618,1.618); // mode 1.0, sd 1.0
}
generated quantities { 
  vector[n_x] beta;
  real beta_0;
  real sigma;

  // Transform to original scale:
  beta = (z_beta ./ x_sd) * y_sd;
  beta_0 = z_beta_0 * y_sd + y_m - sum((z_beta .* x_m) ./ x_sd) * y_sd;
  sigma = z_sigma * y_sd;
}
