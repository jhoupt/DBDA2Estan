data {
  int<lower=1> nCond;
  int<lower=1> nSubj;
  int<lower=1> CondOfSubj[nSubj];
  int<lower=1> nTrlOfSubj[nSubj];

}
parameters {
  vector[nSubj]<lower=0,upper=1> theta;
  real<lower=0,upper=1> omega0;
  vector[nCond]<lower=0,upper=1> omega;
  vector[nCond]<lower=0> kappaMinusTwo;
  real<lower=0,upper=1> modelProb1;
}
transformed parameters {
  vector[nCond]<lower=0> kappa;
    kappa = kappaMinusTwo + 2;
  }
}
model {
  real aP;
  real bP;
  real aBeta;
  real bBeta;
  real logProb1;
  real logProb2;

  // Constants for prior 
  aP = 1;
  bP = 1;

  nCorrOfSubj ~ binomial(nTrlOfSubj,theta);
  for ( s in 1:nSubj ) {
    // Use omega[j] for model index 1, omega0 for model index 2:
    aBeta = omega[CondOfSubj[s]] * kappaMinusTwo[CondOfSubj[s]] +1;
    bBeta = (1-omega[CondOfSubj[s]]) * kappaMinusTwo[CondOfSubj[s]] +1;

    logProb1 = log(modelProb1)   + beta_lpdf( theta[s], aBeta, bBeta );

    aBeta =    omega0  * kappaMinusTwo[CondOfSubj[s]] +1;
    bBeta = (1-omega0) * kappaMinusTwo[CondOfSubj[s]] +1;

    logProb2 = log1m(modelProb1) + beta_lpdf( theta[s], aBeta, bBeta );

    //target +=  log_sum_exp(logProb1, logProb2) ; 
    target +=  log(prob1 + prob2) ; 
  }
  for ( j in 1:nCond ) {
    kappaMinusTwo[j] ~ gamma( 2.618 , 0.0809 ); // mode 20 , sd 20
  }

  omega0 ~ beta( aP, bP );
  for ( j in 1:nCond ) {
    omega[j] ~ beta( aP , bP );
  }

  modelProb1 ~ beta(1,1); 
}
