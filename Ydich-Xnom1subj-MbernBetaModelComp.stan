data {
  int<lower=0> n_total;
  int<lower=0,upper=1> y[n_total];
}
transformed data{
  vector<lower=0,upper=1>[2] omega;
  real<lower=0> kappa;

  omega[1] = .25;
  omega[2] = .75;
  kappa = 12;
}
parameters {
  vector<lower=0,upper=1>[2] theta;
  real<lower=0,upper=1> m_prob;
}
transformed parameters{
  vector<lower=0>[2] alpha;
  vector<lower=0>[2] beta;
  alpha = omega * (kappa - 2) + 1;
  beta  = (1 - omega) * (kappa - 2) + 1;
}
model {
  // Uncomment and change parameters to set prior to something other than 
  // uniform.  
  //m_prob ~ beta(1,1);  

  theta ~ beta(alpha,beta);
  for (i in 1:n_total) {
    target +=  log_mix(m_prob, bernoulli_lpmf(y[i] | theta[1]), 
                               bernoulli_lpmf(y[i] | theta[2]));
  }
}
