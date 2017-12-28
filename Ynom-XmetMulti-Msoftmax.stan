data {
  int<lower=0> n_total;
  int<lower=2> n_out;
  int<lower=2> n_x;
  matrix [n_total, n_x] x;
  int<lower=1,upper=n_out> y[n_total];
}
transformed data { 
  // Standardize the data
  matrix[n_total,n_x] zx;
  vector[n_x] x_m;
  vector[n_x] x_sd;

  for (j in 1:n_x) { 
    x_m[j] = mean(col(x,j));
    x_sd[j] = sd(col(x,j));
    zx[,j] = (col(x,j) - x_m[j]) / x_sd[j]; 
  }
}
parameters { 
  matrix[n_out, n_x] zbeta;
}
model { 
  // Priors vague on standardized scale
  for (j in 1:n_x) { 
    zbeta[j] ~ normal(0, 5);
  }

  for (n in 1:n_total) {
    y[n] ~ categorical(softmax(zbeta * x[n]'));
  }
}
generated quantities {
  // Transform zbeta to original scale:
  matrix[n_out, n_x] beta;

  for ( n in 1:n_out ) {
    beta[n,] = zbeta[n,] ./ x_sd';
  }
}
