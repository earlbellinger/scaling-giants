#### Collect scaling relation exponents and performance into a table 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

scaling_dir <- 'scaling_res'
cv_dir <- 'crossval_res'
comp_dir <- 'crossval_comp'

cutoff.1 <- 0.4
cutoff.2 <- 0.9

options(scipen=10000)

formatter <- function(val, unc) {
    unc <- toString(signif(unc, 2))
    if (nchar(unc) == 1) {
        unc <- paste0(unc, '.0')
    }
    if (as.numeric(unc) < 10) {
        if (substr(unc, nchar(unc)-1, nchar(unc)-1) == '0' || # second char 0
            substr(unc, nchar(unc)-1, nchar(unc)-1) == '.' && # second char .
            substr(unc, nchar(unc)-2, nchar(unc)-2) == '0')   # third  char 0
                unc <- paste0(unc, '0')
        dec_places <- nchar(strsplit(toString(unc), '\\.')[[1]][2])
        c(sprintf(val, fmt=paste0('%.', dec_places, 'f')), unc)
    } else {
        val <- sprintf(val, fmt='%.0f')
        if (nchar(unc) > 2) {
            valx <- unlist(strsplit(val, ''))
            valx[length(valx):(length(valx) - nchar(unc) + 3)] <- '0'
            val <- paste0(valx, collapse='')
        }
        c(val, unc)
    }
}

collect_DF <- function(directory, X, ordering=NULL) {
    filenames <- list.files(directory, full.names=T)
    filenames <- filenames[grepl(X, filenames)]
    DF <- do.call(rbind, Map(function(filename) {
            #print(filename)
            read.table(filename, header=1)
        }, filename=filenames))
    if (!is.null(ordering)) DF <- DF[order(DF[[ordering]]),]
    rownames(DF) <- NULL 
    DF
}

for (X in c('age', 'M', 'R')) {
    
    DF   <- collect_DF(scaling_dir, X, 'combo')
    comp <- collect_DF(comp_dir, X)
    cv   <- collect_DF(cv_dir, X)
    
    cv.med <- as.data.frame(do.call(rbind, Map(function(combo) 
            sapply(cv[cv$combo==combo,], median), 
        combo=DF$combo)))
    
    comp.wm <- as.data.frame(do.call(rbind, Map(function(combo) {
            comp. <- comp[comp$combo == combo,]
            x  <- comp.$ratio
            wt <- 1/comp.$e_ratio**2
            xm <- weighted.mean(x, wt, na.rm=T)
            v  <- sum(wt * (x - xm)**2) / sum(wt) 
            data.frame(mean=xm, sd=sqrt(v))
        }, combo=DF$combo)))
    
    DF$mean  <- comp.wm$mean
    DF$sd    <- comp.wm$sd
    DF$r2    <- cv.med$r2
    DF$r2adj <- cv.med$r2adj
    
    if (X == 'age') {
        age.order <- order(-(DF$r2adj * apply(DF[,2:5] != 0, 1, sum)))
    }
    DF <- DF[age.order,]
    
    DF <- DF[DF$r2adj > cutoff.1,]
    
    tabamp <- ' \t & \t '
    outstr <- ''
    for (ii in 1:nrow(DF)) {
        combo <- DF[ii,]$combo 
        
        sys.unc <- signif(DF$sigma_sys[ii], 2)
        ran.unc <- signif(DF$sigma_ran[ii], 2)
        hline <- ' \\\\'
        num_nonzero <- sum(as.logical(DF[ii,][1:5+1]))
        if (ii != nrow(DF) && 
            num_nonzero != sum(as.logical(DF[ii+1,][1:5+1]))) 
                hline <- paste0(hline, '\\hline')
        
        P <- DF[ii,][1:4+1]
        if (any(sapply(P, function(x) x!=0 && abs(x)<0.1))) next 
        Ps <- paste(sapply(P, function(x) 
                ifelse(x==0, '{--}', signif(x, 4))), 
            collapse=tabamp)
        
        meansd <- formatter(DF[ii,]$mean, DF[ii,]$sd)
        meansd <- paste0(meansd[[1]], '(', 
            substr(meansd[[2]], nchar(meansd[[2]])-1, nchar(meansd[[2]])),
            ')')
        r2 <- paste0(signif(DF[ii,]$r2adj, 2))
        while (nchar(r2) < 4) r2 <- paste0(r2, '0')
        if (r2 == '1000') r2 <- '0.99'
        stats <- paste(sys.unc, meansd, r2, sep=tabamp)
        
        rowcolor <- ''
        if (DF[ii,]$r2 >= cutoff.2) rowcolor <- '\\rowcolor{aliceblue} '        
        
        outstr <- paste(outstr, 
                paste0(rowcolor, 
                    paste(combo, Ps, stats, sep=tabamp), hline), 
            sep='\n')
    }
    outstr <- paste0(outstr, '\n')
    cat(outstr)
    cat('\n\n\n')
}
