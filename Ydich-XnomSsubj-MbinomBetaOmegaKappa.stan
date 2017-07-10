data {
  int<lower=1> n_subj;
  int<lower=0> z[n_subj];
  int<lower=0> N[n_subj];
}
parameters {
  vector<lower=0,upper=1>[n_subj] theta; // individual prob correct
  real<lower=0,upper=1> omega;        // group mode
  real<lower=0> kappa_minus_two;        // group concentration minus two
}
transformed parameters {
  real<lower=0> kappa;  
  kappa = kappa_minus_two + 2;
}
model {
  // Uncomment and change parameters to set prior to something other than 
  // uniform.
  //omega ~ beta(1, 1);

  kappa_minus_two ~ gamma(0.01, 0.01); // mean=1, sd=10 (generic vague)
  // kappa_minus_two ~ gamma(1.105125, 0.1051249);  # mode=1, sd=10 

  theta ~ beta(omega * (kappa - 2) + 1, (1 - omega) * (kappa - 2) + 1);
  z ~ binomial(N, theta);
}
