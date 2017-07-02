
# Stan-Ydich-Xnom1subj-MbernBetaModelComp.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
# Adapted for Stan by Joe Houpt
source("DBDA2E-utilities.R")
fileNameRoot="Stan-Ydich-Xnom1subj-MbernBetaModelComp-" # for output filenames


#===============================================================================

genMCMC = function( data , numSavedSteps=50000 , saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  if ( class(data)=="data.frame" ) {  # If data is a data.frame
    y = myData$y                      # then pull out the column named y
  } else {                            # else
    y = data                          # rename the data as y.
  }
  # Do some checking that data make sense:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=0> Ntotal ;
    int<lower=0,upper=1> y[Ntotal] ;
  }
  transformed data{
    real<lower=0,upper=1> omega[2] ;
    real<lower=0> kappa;
    simplex[2] mPriorProb;

    mPriorProb[1] <- .5;
    mPriorProb[2] <- .5;

    omega[1] <- .25;
    omega[2] <- .75;

    kappa <- 12;
  }
  parameters {
    real<lower=0,upper=1> theta[2] ;
    simplex[2] mProb;
  }
  transformed parameters{
     real<lower=0> alpha[2];
     real<lower=0> beta[2];

     for (m in 1:2) {
       alpha[m] <- omega[m]*(kappa-2)+1 ;
       beta[m]  <- (1-omega[m])*(kappa-2)+1 ;
     }
  }
  model {
    mProb ~ dirichlet(mPriorProb);

    for (m in 1:2) {
       theta[m] ~ beta(alpha[m],beta[m]) ;
    }

    for (i in 1:Ntotal) {
      increment_log_prob( log_sum_exp( log(mProb[1]) + bernoulli_log(y[i], theta[1]), log(mProb[2]) + bernoulli_log(y[i], theta[2]) ) );
    }

  }
  " # close quote for modelString
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  # Option 1: Use single initial value for all chains:
  #  thetaInit = sum(y)/length(y)
  #  initsList = list( theta=thetaInit )
  # Option 2: Use function that generates random values for each chain:
  initsList = function() {
    resampledY = sample( y , replace=TRUE )
    thetaInit = sum(resampledY)/length(resampledY)
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
    resampledY = sample( y , replace=TRUE )
    thetaInit2 = sum(resampledY)/length(resampledY)
    thetaInit2 = 0.001+0.998*thetaInit # keep away from 0,1
    return( list( theta=c(thetaInit, thetaInit2) ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta", "mProb")     # The parameters to be monitored
  burnInSteps = 500            # Stan defaults to iter/2 for overdispersed inits
  nChains = 4                  # nChains should be 2 or more for diagnostics 
  thinSteps = 4                # In Stan there is autocorrelation, so thin

#  stanCpp <- stanc( model_code = modelString ) # Translate to C++
#  stanDso <- stan_model( stanc_ret = stanCpp ) # Compile Stan DSO
  
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( model_code=modelString ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       pars = parameters , # optional
                       chains = nChains ,
                       iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , 
                       warmup = burnInSteps , 
                       thin = thinSteps ,
                       init = initsList ) # optional  
  # Or, accomplish above in one "stan" command; note stanDso is not separate.
  
  # For consistency with JAGS-oriented functions in DBDA2E collection, 
  # convert stan format to coda format:
  codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                               function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
    save( stanFit , file=paste(saveName,"StanFit.Rdata",sep="") )
    save( stanDso , file=paste(saveName,"StanDso.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#------------------------------------------------------------------------------
# THE DATA.

N=9
z=6
y = c( rep(0,N-z) , rep(1,z) )

#------
# Run the model

codaSamples <- genMCMC(y)

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
pM1 = mcmcMat[,"mProb[1]"]
pM2 = mcmcMat[,"mProb[2]"]
thetaM1 = mcmcMat[, "theta[1]"]
thetaM2 = mcmcMat[, "theta[2]"]

# Plot histograms of sampled theta values for each model,
# with pM displayed.
openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,1,2,3),nrow=2,byrow=FALSE) , widths=c(1,2) )
plotPost( pM1 , breaks=seq(0,1,0.05) , cenTend="mean" , xlab="m" , main="Model Index" )
plotPost( thetaM1 , 
          main=bquote( theta*" when m=1" * " ; p(m=1|D)" == .(signif(pM1,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
plotPost( thetaM2 , 
          main=bquote( theta*" when m=2" * " ; p(m=2|D)" == .(signif(pM2,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
saveGraph( file=paste0(fileNameRoot,"Post") , type="eps" )


