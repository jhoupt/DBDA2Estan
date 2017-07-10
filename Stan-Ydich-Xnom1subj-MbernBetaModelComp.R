
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
  if ( class(data)=="list" ) {  # If data is a data.frame
    y = data$y                      # then pull out the column named y
  } else {                            # else
    y = data                          # rename the data as y.
  }
  # Do some checking that data make sense:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    n_total = Ntotal 
  )
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
  parameters = c( "theta", "m_prob")     # The parameters to be monitored
  burnInSteps = 500            # Stan defaults to iter/2 for overdispersed inits
  nChains = 4                  # nChains should be 2 or more for diagnostics 
  thinSteps = 4                # In Stan there is autocorrelation, so thin

#  stanCpp <- stanc( model_code = modelString ) # Translate to C++
#  stanDso <- stan_model( stanc_ret = stanCpp ) # Compile Stan DSO
  
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model(file="Ydich-Xnom1subj-MbernBetaModelComp.stan") 
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
