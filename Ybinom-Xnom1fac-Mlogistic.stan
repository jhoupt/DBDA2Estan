data {
 int<lower=1> Ntotal;
 int<lower=1> NxLvl;
 int<lower=1> N[Ntotal];
 int<lower=0> y[Ntotal];
 int<lower=0> x[Ntotal];
}
parameters {
  real<lower=0,upper=1> mu[Ntotal];
  vector[NxLvl] a;
  real a0;
  real<lower=0> aSigma;
  real<lower=0> kappaMinusTwo;
}
transformed parameters { 
  vector<lower=0,upper=1>[NxLvl] omega;
  real<lower=0> kappa;
  vector[NxLvl] m;
  real b0;
  vector[NxLvl] b;

  omega = inv_logit( a0 + a );
  kappa = kappaMinusTwo + 2;

  // Convert a0,a[] to sum-to-zero b0,b[] :
  m = a0 + a; // cell means 
  b0 = mean( m );
  b = m - b0;
}
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ binomial( N[i], mu[i] );
    mu[i] ~ beta( omega[x[i]]*(kappa-2)+1, (1-omega[x[i]])*(kappa-2)+1 );
  }
  a ~ normal( 0.0, aSigma^2 );
  a0 ~ normal( 0.0 , 2^2 );
  aSigma ~ gamma( 1.64 , 0.32 );  // mode=2, sd=4
  kappaMinusTwo ~ gamma( 0.01 , 0.01 );  // mean=1 , sd=10 (generic vague)
}
