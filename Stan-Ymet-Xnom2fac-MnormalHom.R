
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
# Adapted for Stan by Joe Houpt

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x1 = x1 ,
    x2 = x2 ,
    Ntotal = Ntotal ,
    Nx1Lvl = Nx1Lvl ,
    Nx2Lvl = Nx2Lvl ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    agammaShRa = agammaShRa 
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  data { 
    int<lower=0> Ntotal;
    int<lower=1> Nlevels;
    int<lower=1> Nsubjects;
    real y[Ntotal];
    int<lower=1,upper=Nlevels> x[Ntotal];
    int<lower=1,upper=Nsubjects> subj[Ntotal];
  }
  parameters {
    vector[Nlevels] eox;
    vector[Nsubjects] eos;
    vector<lower=0>[Nsubjects] varOfSubj;
    real<lower=0> grandMeanRT;
  }
  transformed parameters {
    vector[Nlevels] effectOfX;
    vector[Nsubjects] effectOfSubj;
    vector<lower=0>[Ntotal] muTrial; 
    vector<lower=0>[Ntotal] varTrial; 
    vector<lower=0>[Ntotal] shapeTrial; 
    vector<lower=0>[Ntotal] rateTrial; 

    effectOfSubj <- eos - mean(eos);
    effectOfX    <- eox - mean(eox);

    for ( trial in 1:Ntotal ) {
       muTrial[trial] <- grandMeanRT + effectOfX[ x[trial] ] + effectOfSubj[ subj[trial] ] ; 
       varTrial[trial]<-                                          varOfSubj[ subj[trial] ] ; 
       rate[trial] <- muTrial[trial] / varTrial[trial];
       shape[trial] <- muTrial[trial] * rate[trial];
    }
  }
  model {
    for ( trial in 1:Ntotal ) {
      y[trial] ~ gamma( shape[trial], rate[trial] );
    }

    grandMeanRT ~ gamma(1, .001); 
    eox ~ normal(0,100);
    eos ~ normal(0,100);
    varOfSubj ~ gamma(2, 0.00005);
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rjags)
  require(runjags)
  parameters = c( "b0" ,  "b1" ,  "b2" ,  "b1b2" , "m" , "ySigma" , 
                  "a1SD" , "a2SD" , "a1a2SD" )
  adaptSteps = 1000 
  burnInSteps = 2000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================

smryMCMC = function(  codaSamples , 
                      datFrm=NULL , x1Name=NULL , x2Name=NULL ,
                      x1contrasts=NULL , x2contrasts=NULL , x1x2contrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
    x1levels = levels(as.factor(datFrm[,x1Name]))
    x2levels = levels(as.factor(datFrm[,x2Name]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "b1b2[12,34]" then pull out b1b2, 12 and 34:
      strparts = unlist( strsplit( parName , "\\[|,|\\]"  ) )
      # if there are only the param name and a single index:
      if ( length(strparts)==2 ) { 
        # if param name refers to factor 1:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="1" ) { 
          thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])] )
        }
        # if param name refers to factor 2:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="2" ) { 
          thisRowName = paste( thisRowName , x2levels[as.numeric(strparts[2])] )
        }
      }
      # if there are the param name and two indices:
      if ( length(strparts)==3 ) { 
        thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])], 
                             x2levels[as.numeric(strparts[3])] )
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # All contrasts:
  if ( !is.null(x1contrasts) | !is.null(x2contrasts) | !is.null(x1x2contrasts) ) {
    if ( is.null(datFrm) | is.null(x1Name) | is.null(x2Name) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      x1 = as.numeric(as.factor(datFrm[,x1Name]))
      x1levels = levels(as.factor(datFrm[,x1Name]))
      x2 = as.numeric(as.factor(datFrm[,x2Name]))
      x2levels = levels(as.factor(datFrm[,x2Name]))
      # x1 contrasts:
      if ( !is.null(x1contrasts) ) {
        for ( cIdx in 1:length(x1contrasts) ) {
          thisContrast = x1contrasts[[cIdx]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x1levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x1levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b1[",1:length(x1levels),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      }
      # x2 contrasts:
      if ( !is.null(x2contrasts) ) {
        for ( cIdx in 1:length(x2contrasts) ) {
          thisContrast = x2contrasts[[cIdx]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x2levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x2levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b2[",1:length(x2levels),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      }
      # interaction contrasts:
      if ( !is.null(x1x2contrasts) ) {
        for ( cIdx in 1:length(x1x2contrasts) ) {
          thisContrast = x1x2contrasts[[cIdx]]
          # factor A contrast:
          facAcontrast = thisContrast[[1]]
          facAcontrastLeft = facAcontrast[[1]]
          facAcontrastRight = facAcontrast[[2]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( facAcontrastLeft ) ) { 
            left = left | x1levels==facAcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facAcontrastRight ) ) { 
            right = right | x1levels==facAcontrastRight[nIdx]
          }
          right = normalize(right)
          facAcontrastCoef = left-right
          # factor B contrast:
          facBcontrast = thisContrast[[2]]
          facBcontrastLeft = facBcontrast[[1]]
          facBcontrastRight = facBcontrast[[2]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( facBcontrastLeft ) ) { 
            left = left | x2levels==facBcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facBcontrastRight ) ) { 
            right = right | x2levels==facBcontrastRight[nIdx]
          }
          right = normalize(right)
          facBcontrastCoef = left-right
          # interaction contrast:
          contrastCoef = cbind( as.vector( 
            outer( facAcontrastCoef , facBcontrastCoef ) ) )
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,5)=="b1b2["] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste( paste(thisContrast[[1]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[1]][[2]],collapse=""),sep=""),
                   ".x.",
                   paste( paste(thisContrast[[2]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[2]][[2]],collapse=""),sep=""),
                   sep="") )
        }
      }
      
    }
    
  }
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================


plotMCMC = function( codaSamples , 
                     datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                     x1contrasts=NULL , 
                     x2contrasts=NULL , 
                     x1x2contrasts=NULL ,
                     saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  # Display data with posterior predictive distributions
  for ( x2idx in 1:length(x2levels) ) {
    openGraph(width=1.2*length(x1levels),height=5)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    plot(-10,-10, 
         xlim=c(0.2,length(x1levels)+0.1) , 
         xlab=paste(x1Name,x2Name,sep="\n") , 
         xaxt="n" , ylab=yName ,
         ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(x1levels) , tick=FALSE ,
          lab=paste( x1levels , x2levels[x2idx] , sep="\n" ) )
    for ( x1idx in 1:length(x1levels) ) {
      xPlotVal = x1idx #+ (x2idx-1)*length(x1levels)
      yVals = y[ x1==x1idx & x2==x2idx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) , 
              yVals , pch=1 , cex=1.5 , col="red" )
      chainSub = round(seq(1,chainLength,length=20))
      for ( chnIdx in chainSub ) {
        m = mcmcMat[chnIdx,paste("m[",x1idx,",",x2idx,"]",sep="")]
        s = mcmcMat[chnIdx,"ySigma"]
        normlim = qnorm( c(0.025,0.975) )
        yl = m+normlim[1]*s
        yh = m+normlim[2]*s
        ycomb=seq(yl,yh,length=201)
        ynorm = dnorm(ycomb,mean=m,sd=s)
        ynorm = 0.67*ynorm/max(ynorm)
        lines( xPlotVal-ynorm , ycomb , col="skyblue" ) 
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",x2levels[x2idx]), type=saveType)
    }
  } # end for x2idx
  
  # plot contrasts:
  if ( !is.null(x1contrasts) | !is.null(x2contrasts) | !is.null(x1x2contrasts) ) {
    if ( is.null(datFrm) | is.null(x1Name) | is.null(x2Name) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      x1 = as.numeric(as.factor(datFrm[,x1Name]))
      x1levels = levels(as.factor(datFrm[,x1Name]))
      x2 = as.numeric(as.factor(datFrm[,x2Name]))
      x2levels = levels(as.factor(datFrm[,x2Name]))
      # x1 contrasts:
      if ( !is.null(x1contrasts) ) {
        for ( cIdx in 1:length(x1contrasts) ) {
          thisContrast = x1contrasts[[cIdx]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x1levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x1levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b1[",1:length(x1levels),"]",sep="")] 
                           %*% contrastCoef )
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference" ,
                               main=paste0( 
                                 paste(thisContrast[[1]],collapse="."), 
                                 "\nvs\n",
                                 paste(thisContrast[[2]],collapse=".") ),
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          saveGraph(file=paste0(fileNameRoot,
                                paste(thisContrast[[1]],collapse="."), 
                                "-vs-",
                                paste(thisContrast[[2]],collapse=".")),
                    type=graphFileType)
        }
      }
      # x2 contrasts:
      if ( !is.null(x2contrasts) ) {
        for ( cIdx in 1:length(x2contrasts) ) {
          thisContrast = x2contrasts[[cIdx]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x2levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x2levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b2[",1:length(x2levels),"]",sep="")] 
                           %*% contrastCoef )
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference" ,
                               main=paste0( 
                                 paste(thisContrast[[1]],collapse="."), 
                                 "\nvs\n",
                                 paste(thisContrast[[2]],collapse=".") ),
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          saveGraph(file=paste0(fileNameRoot,
                                paste(thisContrast[[1]],collapse="."), 
                                "-vs-",
                                paste(thisContrast[[2]],collapse=".")),
                    type=graphFileType)
        }
      }
      # interaction contrasts:
      if ( !is.null(x1x2contrasts) ) {
        for ( cIdx in 1:length(x1x2contrasts) ) {
          thisContrast = x1x2contrasts[[cIdx]]
          # factor A contrast:
          facAcontrast = thisContrast[[1]]
          facAcontrastLeft = facAcontrast[[1]]
          facAcontrastRight = facAcontrast[[2]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( facAcontrastLeft ) ) { 
            left = left | x1levels==facAcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facAcontrastRight ) ) { 
            right = right | x1levels==facAcontrastRight[nIdx]
          }
          right = normalize(right)
          facAcontrastCoef = left-right
          # factor B contrast:
          facBcontrast = thisContrast[[2]]
          facBcontrastLeft = facBcontrast[[1]]
          facBcontrastRight = facBcontrast[[2]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( facBcontrastLeft ) ) { 
            left = left | x2levels==facBcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facBcontrastRight ) ) { 
            right = right | x2levels==facBcontrastRight[nIdx]
          }
          right = normalize(right)
          facBcontrastCoef = left-right
          # interaction contrast:
          contrastCoef = cbind( as.vector( 
            outer( facAcontrastCoef , facBcontrastCoef ) ) )
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,5)=="b1b2["] 
                           %*% contrastCoef )
          mainName = paste0(
            paste0(paste(thisContrast[[1]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[1]][[2]],collapse=".")),
            "\n(x)\n",
            paste0(paste(thisContrast[[2]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[2]][[2]],collapse=".")))
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference of Differences" ,
                               main=mainName ,
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          fileNameSuffix = paste0(
            paste0(paste(thisContrast[[1]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[1]][[2]],collapse=".")),
            "(x)",
            paste0(paste(thisContrast[[2]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[2]][[2]],collapse=".")))
          saveGraph(file=paste0(fileNameRoot,fileNameSuffix),
                    type=graphFileType)
        }
      } # end of interaction contrasts
    }  
  } # end plot contrasts
}

# 
# #==============================================================================
# # Do NHST ANOVA:
# 
# theData = data.frame( y=y , x1=factor(x1,labels=x1levels) ,
#                             x2=factor(x2,labels=x2levels) )
# openGraph(width=7,height=7)
# interaction.plot( theData$x1 , theData$x2 , theData$y , type="b" )
# #saveGraph( file=paste(fileNameRoot,"DataPlot",sep="") , type="eps" )
# aovresult = aov( y ~ x1 * x2 , data = theData )
# cat("\n------------------------------------------------------------------\n\n")
# print( summary( aovresult ) )
# cat("\n------------------------------------------------------------------\n\n")
# print( model.tables( aovresult , type = "effects", se = TRUE ) , digits=3 )
# cat("\n------------------------------------------------------------------\n\n")
# 
# #==============================================================================
