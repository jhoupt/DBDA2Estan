# Stan-Ynom-XmetMulti-Msoftmax.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
# Adapted for Stan by Joe Houpt

source("DBDA2E-utilities-stan.R")

#===========================================================================

genMCMC = function(data, xName, yName, 
                   numSavedSteps=10000, thinSteps=1, saveName=NULL,
                   nChains=nChainsDefault) { 
  require(rstan)
  #-------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = as.numeric(y) ,
    n_x = dim(x)[2] ,
    n_total = dim(x)[1] ,
    n_out = length( unique( y ) ) # should be same as max(y)
  )
  #-------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let Stan do it...
  #-------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "zbeta" , "beta" )
  adaptSteps = 500  # Number of steps to "tune" the samplers
  burnInSteps = 1000

  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( file="Ynom-XmetMulti-Msoftmax.stan" ) 
  stanFit <- sampling(object = stanDso,
                      data = dataList,
                      pars = parameters,
                      chains = nChains,
                      warmup = burnInSteps,
                      iter = (ceiling(numSavedSteps/nChains)*thinSteps + 
                              burnInSteps)
                      )
  codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                   function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function
#===========================================================================
