#### Plot cross-validated scaling relations 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 

set.seed(42)
crossval_dir <- 'crossval_comp'

variables <- expression(nu['max'], Delta*nu, T['eff'], '[Fe/H]')
combo_sets <- list(c(1, 2, 3, 4), 
    c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4), 
    c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4), 
    c(1), c(2), c(3), c(4))

crossval_filenames <- list.files(crossval_dir)
age_filenames <- crossval_filenames[grepl('age', crossval_filenames)]
M_filenames   <- crossval_filenames[grepl('M',   crossval_filenames)]
R_filenames   <- crossval_filenames[grepl('R',   crossval_filenames)]

age.DF <- do.call(rbind, Map(function(filename)
        read.table(file.path(crossval_dir, filename), header=1), 
    filename=age_filenames))
age.DF <- age.DF[with(age.DF, order(combo, rem_idx)),]

M.DF <- do.call(rbind, Map(function(filename)
        read.table(file.path(crossval_dir, filename), header=1), 
    filename=M_filenames))
M.DF <- M.DF[with(M.DF, order(combo, rem_idx)),]

R.DF <- do.call(rbind, Map(function(filename)
        read.table(file.path(crossval_dir, filename), header=1), 
    filename=R_filenames))
R.DF <- R.DF[with(R.DF, order(combo, rem_idx)),]

cex = 0.66
height = 4.17309*1.5*0.74
width = 4.17309*1.3*1.2

plot_scaling <- function(DF, xlab='Age', majorn=2,
        pch=21, label='', tyticks=T, byticks=T,
        xlim=c(0.4, 2.8), resid_lim=c(-0.3, 0.3), solar_val=1, color='darkblue', 
        ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=-0.28) {
    
    x  <- DF$ratio
    wt <- 1/DF$e_ratio**2
    xm <- weighted.mean(x, wt, na.rm=T)
    v  <- sum(wt * (x - xm)**2) / sum(wt) 
    msd <- formatter(xm, sqrt(v))
    DF <- DF[c(T,F,F),]
    set.seed(42)
    DF <- DF[sample.int(nrow(DF)),]
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.2, cex=text.cex)
    
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim, n=majorn)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=1, lwd=par()$lwd, lty=2)
    
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('topleft', bty='n', cex=1.6*text.cex, pch=NA, legend=label,
        inset=c(-0.05, -0.01))
    par(family=font)
    
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[tyticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext(paste("Scaling Age [Gyr]"), 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    ## lower panel
    
    ylim <- resid_lim 
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=0, lwd=par()$lwd, lty=2)
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('bottomleft', bty='n', cex=1.3*text.cex, pch=NA, 
        inset=c(-0.03, -0.01),
        legend=paste("Mean:", msd[1], '±', msd[2]))
    par(family=font)
    
    rect(xlim[1], ylim[1], xlim[2],1.05 *ylim[1], col='white', border=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex,
        family=font, majorn=2, labels=T, lwd.ticks=par()$lwd, lwd=par()$lwd)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[byticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext("Rel. Diff.", 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext(xlab, 1, 2, outer=F, cex=text.cex)
}

library(RColorBrewer)
col.pal <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")
col.pal <- adjustcolor(colorRampPalette(col.pal)(8))
col.pal[4] <- 'forestgreen' 
col.pal[6] <- 1
col.pal[3] <- orange
col.pal[5] <- blue
col.pal <- c(col.pal[1:7], 1, col.pal[8], c(1,1,1,1,1,1,1,1,1))

for (combo in sort(unique(age.DF$combo))) { 
    age.DF. <- age.DF[age.DF$combo == combo,]
    
    label <- variables[combo_sets[[combo+1]]]
    
    make_plots(plot_scaling, paste0('age_', combo), 
        filepath=file.path('plots', 'crossval_comp'), 
        slides=F, make_png=F, wide=F, tall=F,
        paper_pdf_height=height,
        paper_pdf_width=width,
        cex.paper=cex,
        DF=age.DF., 
        label=label,
        xlab='Age from Grid-Based Modeling [Gyr]',
        tyticks=c(F, T, T),
        majorn=3,
        xlim=c(0, 13), 
        resid_lim=c(-1.6, 1.6), 
        solar_val=4.569,
        color=col.pal[combo+1])
}


plot_scaling <- function(DF, xlab='Mass', majorn=2,
        pch=21, label='', tyticks=T, byticks=T,
        xlim=c(0.4, 2.8), resid_lim=c(-0.3, 0.3), solar_val=1, color='darkblue', ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=-0.28) {
    
    x  <- DF$ratio
    wt <- 1/DF$e_ratio**2
    xm <- weighted.mean(x, wt, na.rm=T)
    v  <- sum(wt * (x - xm)**2) / sum(wt) 
    msd <- formatter(xm, sqrt(v))
    DF <- DF[c(T,F,F),]
    set.seed(42)
    DF <- DF[sample.int(nrow(DF)),]
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.2, cex=text.cex)
    
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim, n=majorn)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=1, lwd=par()$lwd, lty=2)
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('topleft', bty='n', cex=1.6*text.cex, pch=NA, legend=label,
        inset=c(-0.05, -0.01))
    par(family=font)
    
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[tyticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext(paste("Scaling Value [solar units]"), 2, 2.5, 
        outer=F, las=0, cex=text.cex)
    
    
    ## lower panel
    
    ylim <- resid_lim 
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=0, lwd=par()$lwd, lty=2)
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('bottomleft', bty='n', cex=1.3*text.cex, pch=NA, 
        inset=c(-0.03, -0.01),
        legend=paste("Mean:", msd[1], '±', msd[2]))
    par(family=font)
    
    rect(xlim[1], ylim[1], xlim[2],1.05 *ylim[1], col='white', border=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex,
        family=font, majorn=2, labels=T, lwd.ticks=par()$lwd, lwd=par()$lwd)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[byticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext("Rel. Diff.", 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext(xlab, 1, 2, outer=F, cex=text.cex)
}


for (combo in sort(unique(M.DF$combo))) { 
    M.DF. <- M.DF[M.DF$combo == combo,]
    
    label <- variables[combo_sets[[combo+1]]]
    
    make_plots(plot_scaling, paste0('M_', combo), 
        filepath=file.path('plots', 'crossval_comp'), 
        slides=F, make_png=F, wide=F, tall=F,
        paper_pdf_height=height,
        paper_pdf_width=width,
        cex.paper=cex,
        DF=M.DF., 
        label=label,
        xlab='Mass from Grid-Based Modeling [solar units]',
        tyticks=c(T, T),
        majorn=3,
        xlim=c(0.9, 2.4), 
        resid_lim=c(-0.8, 0.8), 
        color='darkblue')
    
}


plot_scaling <- function(DF, xlab='Mass', majorn=2,
        pch=21, label='', tyticks=T, byticks=T,
        xlim=c(0.4, 2.8), resid_lim=c(-0.3, 0.3), solar_val=1, color='darkblue', ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=-0.28) {
    
    x  <- DF$ratio
    wt <- 1/DF$e_ratio**2
    xm <- weighted.mean(x, wt, na.rm=T)
    v  <- sum(wt * (x - xm)**2) / sum(wt) 
    msd <- formatter(xm, sqrt(v))
    DF <- DF[c(T,F,F),]
    set.seed(42)
    DF <- DF[sample.int(nrow(DF)),]
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.2, cex=text.cex)
    
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim, n=majorn)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=1, lwd=par()$lwd, lty=2)
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, scaling_X, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('topleft', bty='n', cex=1.6*text.cex, pch=NA, legend=label,
        inset=c(-0.05, -0.01))
    par(family=font)
    
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[tyticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext(paste("Scaling Value [solar units]"), 2, 2.5, 
        outer=F, las=0, cex=text.cex)
    
    ## lower panel
    
    ylim <- resid_lim 
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim, n=majorn)
    
    abline(a=0, b=0, lwd=par()$lwd, lty=2)
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor('white', alpha.f=1), 
        lwd=0.25, col='white'))
    
    with(DF, points(X, ratio, 
        pch=pch, cex=0.9, 
        bg=adjustcolor(color, alpha.f=0.8), 
        lwd=0.25, col='white'))
    
    par(family="Helvetica")
    legend('bottomleft', bty='n', cex=1.3*text.cex, pch=NA, 
        inset=c(-0.03, -0.01),
        legend=paste("Mean:", msd[1], '±', msd[2]))
    par(family=font)
    
    rect(xlim[1], ylim[1], xlim[2],1.05 *ylim[1], col='white', border=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex,
        family=font, majorn=2, labels=T, lwd.ticks=par()$lwd, lwd=par()$lwd)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=majorn, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    axis(2, yticks[byticks], labels=T, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd) 
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    mtext("Rel. Diff.", 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext(xlab, 1, 2, outer=F, cex=text.cex)
}


for (combo in sort(unique(age.DF$combo))) { 
    R.DF. <- R.DF[R.DF$combo == combo,]
    
    label <- variables[combo_sets[[combo+1]]]
    
    make_plots(plot_scaling, paste0('R_', combo), 
        filepath=file.path('plots', 'crossval_comp'), 
        slides=F, make_png=F, wide=F, tall=F,
        paper_pdf_height=height,
        paper_pdf_width=width,
        cex.paper=cex,
        DF=R.DF., 
        label=label,
        xlab='Radius from Grid-Based Modeling [solar units]',
        byticks=c(T, T, T, T, F),
        tyticks=c(T, T, T),
        majorn=3,
        xlim=c(3, 14), 
        resid_lim=c(-0.8, 0.8), 
        color=red)
}
