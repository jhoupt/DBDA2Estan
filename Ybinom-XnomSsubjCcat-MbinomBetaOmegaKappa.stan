data {
  int<lower=1> n_subj;
  int<lower=1> n_cat;
  int<lower=0> z[n_subj];
  int<lower=0> c[n_subj];
  int<lower=0> N[n_subj];
}
parameters {
  vector<lower=0,upper=1>[n_subj] theta; // individual prob correct
  vector<lower=0,upper=1>[n_cat] omega;   // group mode
  vector<lower=0>[n_cat] kappa_minus_two; // group concentration minus two
  real<lower=0> kappa_minus_twoO; 
  real<lower=0,upper=1> omegaO;
}
transformed parameters {
  real<lower=0> kappaO;  
  vector<lower=0>[n_cat] kappa;  

  kappaO = kappa_minus_twoO + 2;
  kappa = kappa_minus_two + 2;
}
model {

  omegaO ~ beta(1.0, 1.0);
  //omegaO ~ beta(1.025, 1.075); // mode=0.25 , concentration=2.1

  kappa_minus_twoO ~ gamma(0.01, 0.01);  // mean=1 , sd=10 (generic vague)
  //kappa_minus_twoO ~ gamma(1.01005, 0.01005012);  // mode=1 , sd=100
  //kappa_minus_twoO ~ gamma(1.105125, 0.1051249);  // mode=1 , sd=10
  //kappa_minus_twoO ~ gamma(1.105125, 0.01051249);  // mode=10 , sd=100

  omega ~ beta(omegaO * (kappaO - 2) + 1, (1 - omegaO) * (kappaO - 2) + 1);
  kappa_minus_two ~ gamma(0.01, 0.01);// mean=1 , sd=10 (generic vague)

  z ~ binomial(N, theta);
  theta ~ beta(omega[c] .* (kappa[c] - 2) + 1, 
               (1 - omega[c]) .* (kappa[c] - 2) + 1); 
}
