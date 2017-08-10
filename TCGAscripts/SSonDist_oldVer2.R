rm(list=ls()); gc();

altType <- "SNVs"; #SNVs, CNAs
load(sprintf("SSonDistInput_%s.rda", altType))
source("~/work/jeanLib.R")
library(Hmisc)

pSSonParms <- list()
minSamples <- 10;

s <- names(subspaces)[1]
for (s in names(subspaces)) {
    cat(sprintf("Now computing p-values for %s space.\n", s))
    subspace <- as.matrix(subspaces[[s]])
    samplesFile <- sprintf("SSonDistSamples_%s_%s.rda", altType, s)

    if ( file.exists(samplesFile) ) {
        cat(sprintf("%s already exists. Will not resample SSon's.\n",
                    samplesFile))
        load(samplesFile);
    } else {
        i <- 1
        samples <-
            sapply(1:length(nPoss), function(i) {
                ## return(rnorm(1e5)) })
                nPos <- nPoss[i]
                cat(sprintf("Mutation frequency = %d samples.\n", nPos))
                ## Create 1e5 = 10 x 1e4 random CNAs with nPos positive samples
                onfr <-
                    as.numeric(sapply(1:10, function(k) {
                        cat(sprintf("Bundle %d\n", k))
                        Mr <- t(sapply(1:1e4, function(j) {
                            return(sample(c(rep(1, nPos), rep(0, ncol(E) - nPos))))
                        }))
                        Xr <- E %*% t( Mr / nPos - (1-Mr) / (ncol(Mr) - nPos) )
                        
                        XrbackProj <- subspace %*% t( t(Xr) %*% subspace )
                        
                        ## resr <- sapply(1:ncol(Xr), function(j) {
                        ##     ##residuals: SS off Pareto front
                        ##     sum((Xr[,j] - XrbackProj[,j])^2)
                        ## })
                        onfr <- sapply(1:ncol(Xr), function(j) {
                            ##SS on Parto front
                            sum(XrbackProj[,j]^2)
                        })
                        ## apply(Xr, 2, function(x) { sum(x^2) }) - resr - onfr
                        return(onfr)
                    }))
                return(onfr)
            })
        save(samples, file=samplesFile)
    }

    ## Fit the reverse ecdf
    i <- 1
    pdf(sprintf("SSonDistRevCumApprox_%s_%s.pdf", altType, s),
        height=3*3, width=5*3)
    par(mfrow=c(3,5))
    lineFits <- sapply(1:length(nPoss), function(i) {
        nPos <- nPoss[i]
        onfr <- samples[,i]
        
        myecdf <- Ecdf(onfr, what="1-F", log="y",
                       main=sprintf("n = %d", nPos),
                       xlab="SSon", ylab="1-F", lwd=2, col="grey");
        abline(h=minSamples / length(onfr), lty=2)
        abline(v=median(onfr), lty=2)

        ## Use only the portion of the exponantial tail where we have
        ## at least minSamples to estimate the p-value
        sel <- myecdf[["y"]]>minSamples / length(onfr) &
            myecdf[["x"]]>median(onfr)
        toFit <- cbind(myecdf[["x"]], log10(myecdf[["y"]]))[sel,]
        lineFit <- getLine(toFit)
        abline(lineFit, col="red", lwd=2)
        
        names(lineFit) <- c("a", "b")
        return(lineFit)
    })
    plot(nPoss, lineFits[1,], xlab="# positive samples",
         ylab="a (offset)", log="x")
    plot(nPoss, lineFits[2,], xlab="# positive samples", ylab="b (slope)",
         log="x")
    dev.off();

    pSSonParms[[s]] <- lineFits;
    save(nPoss, pSSonParms,
         file=sprintf("SSonDistRevCumApprox_%s.rda", altType))
}


