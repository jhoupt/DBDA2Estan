
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!


source("Stan-Ydich-Xnom1subj-MbernBetaModelComp.R")
fileNameRoot="Stan-Ydich-Xnom1subj-MbernBetaModelComp-" # for output filenames

#------------------------------------------------------------------------------
# THE DATA.

N=9
z=6
y = c( rep(0,N-z) , rep(1,z) )
dataList = list(
  y = y ,
  N = N 
)

#------------------------------------------------------------------------------
# Generate the MCMC chain:
codaSamples <- genMCMC(dataList, numSavedSteps=5000, saveName=fileNameRoot)

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
pM1 = mcmcMat[,"m_prob"]
pM2 = 1-pM1
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


