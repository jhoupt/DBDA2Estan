data {
  int<lower=0> n_total;
  int<lower=1> dim_x; 
  real x[n_total,dim_x];
  int<lower=0,upper=1> y[n_total];
}
transformed data { 
  vector[dim_x] xm;
  vector[dim_x] xsd;
  vector[dim_x] zx[n_total];

  // Standardize the values of X
  for (j in 1:dim_x) {
    xm[j] = mean(x[,j]);
    xsd[j] = sd(x[,j]);
    for (i in 1:n_total) {
      zx[i,j] = (x[i,j] - xm[j]) / xsd[j];
    }
  }
}
parameters { 
   real zbeta0;
   vector[dim_x] zbeta;
}
model {
  for ( i in 1:n_total ) {
    // In Stan, inv_logit is logistic:
    y[i] ~ bernoulli(inv_logit(zbeta0 + dot_product(zbeta, zx[i])));
  }
  // Priors vague on standardized scale:
  zbeta0 ~ normal(0, 2);
  zbeta ~ normal(0, 2);
}
generated quantities { 
  vector[dim_x] beta;
  real beta0;
  // Transform to original scale:
  beta = zbeta ./ xsd; 
  beta0 = zbeta0 - sum(zbeta .* xm ./ xsd);
}
