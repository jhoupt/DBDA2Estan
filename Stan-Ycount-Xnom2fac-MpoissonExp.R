
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R and JAGS, 2nd Edition. Academic Press / Elsevier.
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
  # check that y is non-negative integers:
  if ( !( all.equal(y,round(y,0)) & all(y>=0) ) ) {
    stop("** Y values must be counts, i.e., non-negative integers **")
  }
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  Ncell = length(y) # number of rows of data, not sum(y) 
  # For prior distributions:
  yLogMean = log(sum(y)/(Nx1Lvl*Nx2Lvl)) 
  yLogSD = log(sd(c(rep(0,Ncell-1),sum(y)))) 
  agammaShRa = unlist( gammaShRaFromModeSD( mode=yLogSD , sd=2*yLogSD ) )
  showSDprior = TRUE
  if ( showSDprior ) {
    show( agammaShRa )
    openGraph(height=5,width=7)
    xv = seq(0,yLogSD/2+2*2*yLogSD,length=501)
    plot( xv , dgamma( xv , shape=agammaShRa[1] , rate=agammaShRa[2] ) ,
          xlab="SD" , ylab="p(SD)" , main="Prior on SD parameters" , type="l" ,
          col="skyblue" , lwd=3 )
  }
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x1 = x1 ,
    x2 = x2 ,
    Ncell = Ncell ,
    Nx1Lvl = Nx1Lvl ,
    Nx2Lvl = Nx2Lvl ,
    yLogMean = yLogMean ,
    yLogSD = yLogSD ,
    agammaShRa = agammaShRa 
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> Ncell;
    int<lower=0> y[Ncell];
    int x1[Ncell];
    int x2[Ncell];
    int<lower=1>Nx1Lvl;
    int<lower=1>Nx2Lvl;
    real yLogMean;
    real yLogSD;
    real agammaShRa[2];
  }
  parameters {
    real a0;
    real a1[Nx1Lvl];
    real a2[Nx2Lvl];
    real<lower=0> a1SD;
    real<lower=0> a2SD;
    real a1a2[Nx1Lvl,Nx2Lvl];
    real<lower=0> a1a2SD;
  }
  transformed parameters { 
    real<lower=0> lambda[Ncell];
    matrix[Nx1Lvl,Nx2Lvl] m;
    real b0;
    real b1[Nx1Lvl];
    real b2[Nx2Lvl];
    real b1b2[Nx1Lvl,Nx2Lvl];
    matrix[Nx1Lvl,Nx2Lvl] expm;
    matrix[Nx1Lvl,Nx2Lvl] ppx1x2p;
    real<lower=0> ppx1p[Nx1Lvl];
    real<lower=0> ppx2p[Nx2Lvl];

    for ( i in 1:Ncell ) {
      lambda[i] <- exp( a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]] );
    }


    # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2]; # cell means 
    } }

    b0 <- mean( m );
    for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( row(m,j1) ) - b0; }
    for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( col(m,j2) ) - b0; }
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] );
    } }    

    # Compute predicted proportions:
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      expm[j1,j2] <- exp(m[j1,j2]);
    } }
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      ppx1x2p[j1,j2] <- expm[j1,j2]/sum(expm);
    } }
    for ( j1 in 1:Nx1Lvl ) { ppx1p[j1] <- sum(row(ppx1x2p,j1)); }
    for ( j2 in 1:Nx2Lvl ) { ppx2p[j2] <- sum(col(ppx1x2p,j2)); }

  }
  model {
    for ( i in 1:Ncell ) {
      y[i] ~ poisson( lambda[i] );
    }
    a0 ~ normal( yLogMean , (yLogSD*2)^2 ) ;

    for ( j1 in 1:Nx1Lvl ) { a1[j1] ~ normal( 0.0 , a1SD^2 ); }
    a1SD ~ gamma(agammaShRa[1],agammaShRa[2]);

    for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ normal( 0.0 , a2SD^2 ); }
    a2SD ~ gamma(agammaShRa[1],agammaShRa[2]) ;

    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      a1a2[j1,j2] ~ normal( 0.0 , a1a2SD^2 );
    } }
    a1a2SD ~ gamma(agammaShRa[1],agammaShRa[2]) ;
  }
  " # close quote for modelstring
  writeLines(modelString,con="model.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let Stan do it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rstan)
#  parameters = c( "ppx1x2p" , "ppx1p" , "ppx2p" , 
#                 "b0" ,  "b1" ,  "b2" ,  "b1b2" , "m" ,  
#                 "a1SD" , "a2SD" , "a1a2SD" )
  adaptSteps = 1000 
  burnInSteps = 2000 

  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( model_code=modelString ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       #pars = parameters , # optional
                       chains = nChains ,
                       iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , 
                       warmup = burnInSteps , 
                       thin = thinSteps,
                       init = 'random')
                       
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
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,4)=="b1b2"] 
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
  Nx1levels = length(x1levels)
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Nx2levels = length(x2levels)
  # Display data with posterior predictive distributions
  openGraph(width=2.25*Nx2levels,height=1.5*Nx1levels)
  par( mar=c(3.5,2.5,3.0,1.5) , mgp=c(2,0.7,0) )
  layout(matrix(1:(Nx1levels*Nx2levels),nrow=Nx1levels,byrow=TRUE))
  xLim = range(mcmcMat[,grep("^ppx1x2p",colnames(mcmcMat))])
  for ( x1Idx in 1:Nx1levels ) {
    for ( x2Idx in 1:Nx2levels ) {
      cellN = datFrm[ datFrm[,x1Name]==x1levels[x1Idx] 
                      & datFrm[,x2Name]==x2levels[x2Idx] , yName ]
      totalN = sum(y)
      paramCol = paste0("ppx1x2p[",x1Idx,",",x2Idx,"]")
      plotPost( mcmcMat[,paramCol] , xlab="Proportion" , cex.lab=1.0 ,
                main=paste0( x1Name,":",x1levels[x1Idx] ,"  ",
                             x2Name,":",x2levels[x2Idx] ,"\n",
                             "N=",cellN ) , cenTend="mode" ,
                xlim=xLim , border="skyblue" , HDItextPlace=0.9 )
      points(  cellN/totalN ,  0 , pch=17 , col="red" , cex=2 )
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"PostPred"), type=saveType)
  }
  
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
          compInfo = plotPost( postContrast , xlab="Beta Deflect. Diff." ,
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
          compInfo = plotPost( postContrast , xlab="Beta Deflect. Diff." ,
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
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,4)=="b1b2"] 
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
          compInfo = plotPost( postContrast , xlab="Beta Deflect. Diff. of Diff." ,
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
