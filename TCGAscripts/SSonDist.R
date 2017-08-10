rm(list=ls()); gc();
source("~/work/jeanLib.R")
library(Hmisc)

altType <- "SNVs"; #SNVs, CNAs
## altType <- "CNAs"; #SNVs, CNAs
nBundles <- 2; #FIXME 2 for testing
nBundles <- 10; #10 for real
load(sprintf("SSonDistInput_%s.rda", altType))

## maxFreq <- 1e3;
## if ( max(nPoss) > maxFreq ) { maxFreq <- max(nPoss) }

pSSonParms <- list()
minSamples <- 10;
nPoss <- floor(exp(seq(log(minSamples), log(.99 * ncol(E)), len=12)))

##################################################
## Norms don't depend on subspace

samplesFile <- sprintf("SSonDistSamples_%s_%s.rda", altType, "norms")

if ( file.exists(samplesFile) ) {
    cat(sprintf("%s already exists. Will not resample.\n",
                samplesFile))
    load(samplesFile);
} else {
    norms <-
        sapply(1:length(nPoss), function(i) {
            nPos <- nPoss[i]
            cat(sprintf("Sample norms for mutation frequency = %d samples.\n", nPos))
            ## Create 1e5 = nBundles x 1e4 random CNAs with nPos positive samples
            k <- 1
            
            norms <-
                as.numeric(
                    sapply(1:nBundles, function(k) {
                        cat(sprintf("Bundle %d\n", k))
                        Mr <- t(sapply(1:1e4, function(j) {
                            return(sample(c(rep(1, nPos), rep(0, ncol(E) - nPos))))
                        }))
                        Xr <- E %*% t( Mr / nPos - (1-Mr) / (ncol(Mr) - nPos) )
                        
                        return(
                            sapply(1:ncol(Xr), function(j) {
                                sqrt(sum(Xr[,j]^2))
                            }))
                    })
                )
            return(norms)
        })
   save(norms, file=samplesFile);
}

i <- 1
pdf(sprintf("SSonDist_normDistRevCumApprox_%s.pdf", altType),
    height=3*3, width=5*3)
par(mfrow=c(3,5))
i <- 1
SSlineFits <- sapply(1:length(nPoss), function(i) {
    nPos <- nPoss[i]
    myNorms <- norms[,i]
    
    myecdf <- Ecdf(myNorms^2, what="1-F", log="y",
                   main=sprintf("n = %d", nPos),
                   xlab="norm", ylab="1-F", lwd=2, col="grey");
    abline(h=minSamples / length(myNorms), lty=2)
    abline(v=median(myNorms^2), lty=2)
    
    ## Use only the portion of the exponantial tail where we have
    ## at least minSamples to estimate the p-value
    sel <- myecdf[["y"]]>minSamples / length(myNorms) &
        myecdf[["x"]]>median(myNorms^2)
    toFit <- cbind(myecdf[["x"]], log10(myecdf[["y"]]))[sel,]
    lineFit <- getLine(toFit)
    abline(lineFit, col="red", lwd=2)
    
    names(lineFit) <- c("a", "b")
    return(lineFit)
})
plot(nPoss, SSlineFits[1,], xlab="# positive samples",
     ylab="a (offset)", log="x")
plot(nPoss, SSlineFits[2,], xlab="# positive samples", ylab="b (slope)",
     log="x")
dev.off();

pSSonParms[["SS"]] <- SSlineFits;
save(nPoss, pSSonParms,
     file=sprintf("SSonDistRevCumApprox_%s.rda", altType))

##################################################
## Ratios depend on subspace but not on mutation frequency
mynPos <- min(nPoss);
mynPos <- max(nPoss);
mynPos <- round(median(nPoss));

pdf(sprintf("SSonDist_ratioDistRevCumApprox_%s.pdf", altType),
    height=4*1, width=4*1)
par(mfrow=c(1,1))

s <- names(subspaces)[1]
for (s in names(subspaces)) {
    cat(sprintf("Now computing p-values for %s space.\n", s))
    subspace <- as.matrix(subspaces[[s]])
    samplesFile <- sprintf("SSonDistSamples_%s_%s.rda", altType, s)

    if ( file.exists(samplesFile) ) {
        cat(sprintf("%s already exists. Will not resample ratios.\n",
                    samplesFile))
        load(samplesFile);
    } else {
        k <- 1
        ratios <-
            as.numeric(sapply(1:nBundles, function(k) {
                cat(sprintf("Bundle %d\n", k))
                Mr <- t(sapply(1:1e4, function(j) {
                    return(sample(c(rep(1, mynPos), rep(0, ncol(E) - mynPos))))
                }))
                Xr <- E %*% t( Mr / mynPos - (1-Mr) / (ncol(Mr) - mynPos) )

                XrbackProj <- subspace %*% t( t(Xr) %*% subspace )
                        
                onfr <- sapply(1:ncol(Xr), function(j) {
                    ##SS on Parto front
                    sum(XrbackProj[,j]^2)
                })
                
                resr <- sapply(1:ncol(Xr), function(j) {
                    ##residuals: SS off Pareto front
                    sum((Xr[,j] - XrbackProj[,j])^2)
                })

                return(onfr / (resr + onfr))
            }))
        save(ratios, file=samplesFile)
    }

    ## Fit the reverse ecdf
    myecdf <- Ecdf(ratios, what="1-F", log="y", pl=F,
                   main=s, ## ylim=c(1e-7, 1),
                   xlab="SSon", ylab="1-F", lwd=2, col="grey")
    ## abline(h=minSamples / length(ratios), lty=2)

    ## Use only the portion of the exponantial tail where we have
    ## at least minSamples to estimate the p-value
    ## FIXME a line isn't a proper fit for sure... maybe a sine or an
    ## exponential? Or a second order polynom?
    sel <- myecdf[["y"]]>0
    sel <- myecdf[["y"]]>minSamples / length(ratios) & myecdf[["y"]]>0
    sel <- myecdf[["y"]]>minSamples / length(ratios) &
        myecdf[["y"]] < 1 & myecdf[["y"]]>0

    plot(myecdf[["x"]], log10(myecdf[["y"]]), main=s,
         xlim=c(0,1), ylim=c(-8,0),
         xlab="fraction of SS on front", ylab="log10 1-F")
    abline(h=log10(minSamples / length(ratios)), lty=2)
    toFit <- cbind(myecdf[["x"]], log10(myecdf[["y"]]))[sel,]
    points(toFit[,1], toFit[,2], col="grey");
    
    toExport <-
        sapply(runif(1e3) * diff(range(toFit[,1])) +
        min(range(toFit[,1])), function(x) { which.min(abs(x -
        toFit[,1])) })
    ## plot(toFit[toExport,1], toFit[toExport,2], col="grey");
    toFit <- toFit[toExport,]

    optRes <-
        optim(par=rep(0, 6), function(par) {
            sum(abs(
                sapply(toFit[,1], function(x) {
                    sum(
                        par * x^(seq(0, length(par)-1))
                    )
                }) - toFit[,2]
            ))
        })
    ## xs <- seq(min(myecdf[["x"]]), max(myecdf[["x"]]), len=100)
    xs <- seq(0, 1, len=100)
    lines(xs,
          sapply(xs, function(x) {
              sum(
                  optRes$par * x^(seq(0, length(optRes$par)-1))
              )
          }), col="blue")

    ## lineFit <- getLine(toFit)
    ## abline(lineFit, col="red", lwd=2)
    ## names(lineFit) <- c("a", "b")
    ## pSSonParms[[s]] <- lineFit;

    pSSonParms[[s]] <- optRes$par;
}
dev.off()

##################################################
## Norms on the front depend both on subspace (we'll only do front)
## and on mutation frequency

samplesFile <- sprintf("SSonDistSamples_%s_%s.rda", altType, "normOnFront")

if ( file.exists(samplesFile) ) {
    cat(sprintf("%s already exists. Will not resample.\n",
                samplesFile))
    load(samplesFile);
} else {
    s <- "Front"
    subspace <- as.matrix(subspaces[[s]])
    
    norms <-
        sapply(1:length(nPoss), function(i) {
            nPos <- nPoss[i]
            cat(sprintf("Sample norms on front for mutation frequency = %d samples.\n", nPos))
            ## Create 1e5 = nBundles x 1e4 random CNAs with nPos positive samples
            k <- 1
            
            norms <-
                as.numeric(
                    sapply(1:nBundles, function(k) {
                        cat(sprintf("Bundle %d\n", k))
                        Mr <- t(sapply(1:1e4, function(j) {
                            return(sample(c(rep(1, nPos), rep(0, ncol(E) - nPos))))
                        }))
                        Xr <- E %*% t( Mr / nPos - (1-Mr) / (ncol(Mr) - nPos) )
                        XrbackProj <-
                            subspace %*% t( t(Xr) %*% subspace )

                        onfr <- sapply(1:ncol(Xr), function(j) {
                            ##SS on Parto front
                            sum(XrbackProj[,j]^2)
                        })
                        
                        return(onfr)
                    })
                )
            return(norms)
        })
   save(norms, file=samplesFile);
}

i <- 1
pdf(sprintf("SSonDist_normOnFrontDistRevCumApprox_%s.pdf", altType),
    height=3*3, width=5*3)
par(mfrow=c(3,5))
i <- 1
SSlineFits <- sapply(1:length(nPoss), function(i) {
    nPos <- nPoss[i]
    myNorms <- norms[,i] ## myNorms is actually sum of squares
    
    myecdf <- Ecdf(myNorms, what="1-F", log="y",
                   main=sprintf("n = %d", nPos),
                   xlab="norm", ylab="1-F", lwd=2, col="grey");
    abline(h=minSamples / length(myNorms), lty=2)
    abline(v=median(myNorms), lty=2)
    
    ## Use only the portion of the exponantial tail where we have
    ## at least minSamples to estimate the p-value
    sel <- myecdf[["y"]]>minSamples / length(myNorms) &
        myecdf[["x"]]>median(myNorms)
    toFit <- cbind(myecdf[["x"]], log10(myecdf[["y"]]))[sel,]
    lineFit <- getLine(toFit)
    abline(lineFit, col="red", lwd=2)
    
    names(lineFit) <- c("a", "b")
    return(lineFit)
})
plot(nPoss, SSlineFits[1,], xlab="# positive samples",
     ylab="a (offset)", log="x")
plot(nPoss, SSlineFits[2,], xlab="# positive samples", ylab="b (slope)",
     log="x")
dev.off();

pSSonParms[["SSonFront"]] <- SSlineFits;

save(nPoss, pSSonParms,
     file=sprintf("SSonDistRevCumApprox_%s.rda", altType))

pSSonParms
