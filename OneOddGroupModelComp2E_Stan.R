# OneOddGroupModelComp2E.R
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
# Adapted for Stan by Joe Houpt
graphics.off()
rm(list=ls(all=TRUE))
source("DBDA2E-utilities.R")
#require(rjags)
#require(runjags)
fileNameRoot="OneOddGroupModelComp2E-" 

#------------------------------------------------------------------------------
# THE DATA.

# Randomly generated fictitious data.
# For each subject, specify the condition s/he was in,
# the number of trials s/he experienced, and the number correct.
npg = 20  # number of subjects per group
ntrl = 20 # number of trials per subject
CondOfSubj = c( rep(1,npg) , rep(2,npg) , rep(3,npg) , rep(4,npg) )
nTrlOfSubj = rep( ntrl , 4*npg )
set.seed(47405)
condMeans = c(.40,.50,.51,.52)
nCorrOfSubj = c( rbinom(npg,ntrl,condMeans[1]) , rbinom(npg,ntrl,condMeans[2]) ,
                 rbinom(npg,ntrl,condMeans[3]) , rbinom(npg,ntrl,condMeans[4]) )
nCond = length(unique(CondOfSubj))
nSubj = length(CondOfSubj)
# jitter the data to be as close as possible to desired condition means:
for ( cIdx in 1:nCond ) {
  nToAdd = round(condMeans[cIdx]*npg*ntrl)-sum(nCorrOfSubj[CondOfSubj==cIdx])
  if ( nToAdd > 0 ) {
    for ( i in 1:nToAdd ) {
      thisNcorr = ntrl
      while ( thisNcorr == ntrl ) {
        randSubjIdx = sample(which(CondOfSubj==cIdx),size=1)
        thisNcorr = nCorrOfSubj[randSubjIdx]
      }
      nCorrOfSubj[randSubjIdx] = nCorrOfSubj[randSubjIdx]+1
    }
  } 
  if ( nToAdd < 0 ) {
    for ( i in 1:abs(nToAdd) ) {
      thisNcorr = 0
      while ( thisNcorr == 0 ) {
        randSubjIdx = sample(which(CondOfSubj==cIdx),size=1)
        thisNcorr = nCorrOfSubj[randSubjIdx]
      }
      nCorrOfSubj[randSubjIdx] = nCorrOfSubj[randSubjIdx]-1
    }
  }
}


show( aggregate( nCorrOfSubj , by=list(CondOfSubj) , FUN=mean ) / ntrl )

# Package the data:
dataList = list(
  nCond = nCond ,
  nSubj = nSubj ,
  CondOfSubj = CondOfSubj ,
  nTrlOfSubj = nTrlOfSubj ,
  nCorrOfSubj = nCorrOfSubj
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Let Stan do it...

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("omega","kappa","omega0","theta","modelProb1")
adaptSteps = 1000            # Number of steps to "tune" the samplers.
burnInSteps = 5000           # Number of steps to "burn-in" the samplers.
nChains = 3                  # Number of chains to run.
numSavedSteps=12000          # Total number of steps in chains to save.
thinSteps=10                 # Number of steps to "thin" (1=keep every step).

# Translate to C++ and compile to DSO:
stanDso <- stan_model( file="OneOddGroupModelComp.stan" ) 
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     pars = parameters , # optional
                     chains = nChains ,
                     iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                              +burnInSteps ) , 
                     warmup = burnInSteps , 
                     thin = thinSteps ,
                     init = "random" ) # optional  
  # Or, accomplish above in one "stan" command; note stanDso is not separate.

# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format:
codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                             function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
save( codaSamples , file=paste(fileNameRoot,"Mcmc.Rdata",sep="") )
save( stanFit , file=paste(fileNameRoot,"StanFit.Rdata",sep="") )
save( stanDso , file=paste(fileNameRoot,"StanDso.Rdata",sep="") )
 


#------------------------------------------------------------------------------- 
# Display diagnostics of chain:

parameterNames = varnames(codaSamples) # get all parameter names
show(parameterNames)
for ( parName in c("modelProb1","omega[1]","omega0","kappa[1]","theta[1]") ) { 
  diagMCMC( codaSamples , parName=parName ,
            saveName=fileNameRoot , saveType="eps" )
}

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

mcmcMat = as.matrix(codaSamples,chains=TRUE)

xLim=c(0.35,0.75)

# Display the model index
pM1 = mcmcMat[, "modelProb1" ]
pM2 = 1 - pM1
string1 =paste("p( Diff Omega M1 | D )=",round(mean(pM1),3),sep="")
string2 =paste("p( Same Omega M2 | D )=",round(mean(pM2),3),sep="")
openGraph(10,4)
nStepsToPlot = 1000
plot( 1:nStepsToPlot , pM1[1:nStepsToPlot] , type="l" , lwd=2 , ylim=c(0,1),
      xlab="Step in Markov chain" , ylab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") , col="skyblue" )
saveGraph(file=paste0(fileNameRoot,"modelProb1"),type="eps")

# Display the omega0 posterior
omega0sample = mcmcMat[, "omega0" ]
openGraph()
#layout( matrix(1:2,nrow=2) )
plotPost( omega0sample , main="Posterior for Same Omega"  ,
      xlab=expression(omega[0]) , xlim=xLim )
saveGraph(file=paste0(fileNameRoot,"Omega0"),type="eps")

# Display the omega[j] posterior
omega1sample = mcmcMat[, "omega[1]" ]
omega2sample = mcmcMat[, "omega[2]" ]
omega3sample = mcmcMat[, "omega[3]" ]
omega4sample = mcmcMat[, "omega[4]" ]
openGraph(10,5)
layout( matrix(1:4,nrow=1,byrow=T) )
plotPost( omega1sample , main="Posterior for Diff Omega" ,
          xlab=expression(omega[1]) , xlim=xLim )
plotPost( omega2sample , main="Posterior for Diff Omega" ,
          xlab=expression(omega[2]) , xlim=xLim )
plotPost( omega3sample , main="Posterior for Diff Omega" ,
          xlab=expression(omega[3]) , xlim=xLim )
plotPost( omega4sample , main="Posterior for Diff Omega" ,
          xlab=expression(omega[4]) , xlim=xLim )
saveGraph(file=paste0(fileNameRoot,"OmegaCond"),type="eps")


# Display the differences of omega[j]'s
omegaSample = rbind( omega1sample , omega2sample , omega3sample , omega4sample )
openGraph(10,5)
layout( matrix(1:6,nrow=2,ncol=3,byrow=T) )
xmin = -0.25
xmax = 0.25
for ( i in 1:3 ) {
    for ( j in (i+1):4 ) {
        plotPost( omegaSample[i,]-omegaSample[j,] , compVal=0.0 ,
                  xlab=bquote(omega[.(i)]-omega[.(j)]) ,
                  #breaks=unique( c( min(c(xmin,omegaSample[i,]-omegaSample[j,])),
                  #          seq(xmin,xmax,len=20),
                  #          max(c(xmax,omegaSample[i,]-omegaSample[j,])) )) ,
                  main="" , xlim=c(xmin,xmax) )
    }
}
saveGraph(file=paste0(fileNameRoot,"OmegaDiff"),type="eps")
