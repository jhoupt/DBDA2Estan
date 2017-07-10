# Stan-Ymet-XmetSsubj-MrobustHierQuadWt.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan 2nd Edition. Academic Press / Elsevier.
# Adapted for Stan by Joe Houpt

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , sName="s" , wName=NULL ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  s = as.numeric(data[,sName])
  if ( !is.null(wName) ) {
    w = data[,wName]
  } else {
    w = rep(1,length(y))
  }
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s ,
    w = w ,
    n_subj = max(s)  , # should equal length(unique(s))
    n_total = length(y)
  )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  
  # Use lm() to find reasonable coefficients overall, then start all individual
  # units and overall at those values.
  # N.B. THIS DOES NOT ALWAYS WORK AND DOES NOT ALWAYS IMPROVE THE MCMC SAMPLE.
  # IF IT'S A PROBLEM, COMMENT OUT THE inits ARGUMENT IN THE run.jags COMMAND.
  zx = ( x - mean(x) ) / sd(x)
  zxsq = zx^2
  zy = ( y - mean(y) ) / sd(y)
  lmInfo = lm( zy ~ zx + zxsq )
  b0init = lmInfo$coef[1]
  b1init = lmInfo$coef[2]
  sigmaInit = sqrt(mean(lmInfo$res^2))
  nuInit = 10 # arbitrary
  initsList = list(
    z_sigma=sigmaInit  ,
    nu=nuInit ,
    z_beta_0_mu=b0init ,
    z_beta_1_mu=b1init ,
    z_beta_0=rep(b0init,max(s)) ,
    z_beta_1=rep(b1init,max(s))
  )
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta_0", "beta_1" , 
                  "beta_0mu", "beta_1_mu" ,
                  "z_beta_0", "zbeta_1" , 
                  "z_beta_0_mu" , "zbeta_1_mu" , 
                  "sigma" , "nu" , 
                  "z_sigma", "z_beta_0_sigma" , "z_beta_1_sigma"  )
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000 
  nChains = 3 
  
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( file="Stan-Ymet-XmetSsubj-MrobustHier.stan" ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       #pars = parameters , # optional
                       #init = initsList , # optional  
                       chains = nChains ,
                       iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , 
                       warmup = burnInSteps , 
                       thin = thinSteps )
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

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , 
                     xName="x" , yName="y" , sName="s" , wName="w" ,
                     compValBeta0=NULL , ropeBeta0=NULL , 
                     compValBeta1=NULL , ropeBeta1=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = data[,xName]
  s = factor(data[,sName])
  nSubj = length(unique(s)) # should be same as max(s)
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  beta0mu = mcmcMat[,"beta0mu"]
  beta1mu = mcmcMat[,"beta1mu"]
  zbeta0mu = mcmcMat[,"zbeta0mu"]
  zbeta1mu = mcmcMat[,"zbeta1mu"]
  sigma = mcmcMat[,"sigma"]
  nu = mcmcMat[,"nu"]
  log10nu = log10(nu)
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0mu , beta1mu , sigma , log10nu )[plotIdx,] ,
           labels=c( expression(mu[beta*0]) , expression(mu[beta*1]) , 
                     expression(mu[beta*2]) , 
                     expression(sigma) ,  expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  # Set up window and layout:
  nPtToPlot = 1000
  plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
  openGraph(width=8,height=8)
  layout( matrix( 1:9 , nrow=3, byrow=TRUE ) )
  par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( beta0mu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(mu[beta*0]) , main=paste("Intercept, Group Level") )
  histInfo = plotPost( beta1mu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(mu[beta*1]) , main=paste("Slope, Group Level") )
  histInfo = plotPost( zbeta0mu , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(zmu[beta*0]) , main=paste("Intercept, Group Level") )
  histInfo = plotPost( zbeta1mu , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(zmu[beta*1]) , main=paste("Slope, Group Level") )
  #plot( beta1mu[plotIdx] , beta0mu[plotIdx] , 
  #      xlab=bquote(mu[beta*1]) , ylab=bquote(mu[beta*0]) ,
  #      col="skyblue" , cex.lab = 1.75 )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValSigma , ROPE=ropeSigma ,
                       xlab=bquote(sigma) , main=paste("Scale, Subj Level") )
  histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=NULL , ROPE=NULL ,
                       xlab=bquote(log10(nu)) , main=paste("Normality, Subj Level") )
  plot( log10nu[plotIdx] , sigma[plotIdx] , 
        xlab=bquote(log10(nu)) ,ylab=bquote(sigma) , 
        col="skyblue" , cex.lab = 1.75 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  nPanels=25
  nPlots = ceiling(nSubj/nPanels)
  for ( plotIdx in 1:nPlots ) {
    openGraph()
    par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
    layout(matrix(1:nPanels,nrow=5,byrow=TRUE))
    xRang = max(x)-min(x)
    yRang = max(y)-min(y)
    xLimMult = 0.1
    yLimMult = 0.1
    xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
    yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
    #for ( sIdx in unique(ceiling(seq(1,nSubj,length=nPanels))) ) { 
    for ( sIdx in ((plotIdx-1)*nPanels+1):min(nSubj,(plotIdx-1)*nPanels+nPanels)) { 
      thisSrows = (as.numeric(s)==sIdx)
      plot( x[thisSrows] , y[thisSrows] , 
            cex=1.0 , lwd=1 , col="black" , xlim=xLim , ylim=yLim ,
            xlab=xName , ylab=yName , cex.lab=1.0 ,
            main=paste0("Unit: ",levels(s)[sIdx]) , 
            cex.main=1.0  )
      # Superimpose a smattering of believable regression lines:
      nPredCurves=30
      xComb = seq(xLim[1],xLim[2],length=301)
      for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
        b0 = mcmcMat[i,paste0("beta0[",sIdx,"]")]
        b1 = mcmcMat[i,paste0("beta1[",sIdx,"]")]
        lines( xComb , b0+b1*xComb , col="skyblue" )
      }
      points( x[thisSrows] , y[thisSrows] , pch=19 )
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPredSubj",plotIdx), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  openGraph()
  par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
  # Plot data values:
  xRang = max(x)-min(x)
  yRang = max(y)-min(y)
  xLimMult = 0.2
  yLimMult = 0.2
  xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
  yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
  plot( x , y , pch="" , cex=1.0 , col="black" , 
        xlim=xLim , ylim=yLim ,
        xlab=xName , ylab=yName , cex.lab=1.0 ,
        main="All Units" , cex.main=1.0  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=70
  xComb = seq(xLim[1],xLim[2],length=301)
  for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
    b0 = mcmcMat[i,paste0("beta0mu")]
    b1 = mcmcMat[i,paste0("beta1mu")]
    lines( xComb , b0+b1*xComb , col="skyblue" )
  }
  for ( sIdx in 1:nSubj ) {
    thisSrows = (as.numeric(s)==sIdx)
    lines( x[thisSrows] , y[thisSrows] , type="o" , pch=19 ) #, pch=sIdx , col=sIdx )
  }
  #
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPredAll",sep=""), type=saveType)
  }
}

#===============================================================================
