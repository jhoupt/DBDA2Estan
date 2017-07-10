data {
  int<lower=1> n_cond;
  int<lower=1> n_subj;
  int<lower=1> cond_of_subj[n_subj];
  int<lower=1> n_trl_of_subj[n_subj];
  int<lower=0> n_corr_of_subj[n_subj];

}
transformed data {
  real a_p;
  real b_p;
  // Constants for prior 
  a_p = 1;
  b_p = 1;

}
parameters {
  vector<lower=0,upper=1>[n_subj] theta;
  real<lower=0,upper=1> omega0;
  vector<lower=0,upper=1>[n_cond] omega;
  vector<lower=0>[n_cond] kappa_minus_two;
  real<lower=0,upper=1> model_prob_1;
}
transformed parameters {
  vector<lower=0>[n_cond] kappa;
  kappa = kappa_minus_two + 2;
}
model {
  // Uncomment and set parameters to something other than 1, 1 for 
  // non-uniform prior
  //model_prob_1 ~ beta(1,1); 

  n_corr_of_subj ~ binomial(n_trl_of_subj,theta);

  for ( s in 1:n_subj ) {
    real a_beta0;
    real b_beta0;
    real a_beta;
    real b_beta;

    // Use omega[j] for model index 1, omega0 for model index 2:
    a_beta = omega[cond_of_subj[s]] * kappa_minus_two[cond_of_subj[s]] + 1;
    b_beta = (1 - omega[cond_of_subj[s]]) * kappa_minus_two[cond_of_subj[s]]
             + 1;

    a_beta0 =    omega0  * kappa_minus_two[cond_of_subj[s]] + 1;
    b_beta0 = (1-omega0) * kappa_minus_two[cond_of_subj[s]] + 1;


    target +=  log_mix(model_prob_1, beta_lpdf(theta[s] | a_beta, b_beta),
                       beta_lpdf(theta[s] | a_beta0, b_beta0));
  }
  kappa_minus_two ~ gamma(2.618, 0.0809); // mode 20 , sd 20

  omega0 ~ beta(a_p, b_p);
  omega ~ beta(a_p, b_p);
}
