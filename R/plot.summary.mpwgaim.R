## Method plot.summary.mpwgaim
##
##' Plot of an object of class  'summary.mpwgaim'
##'
##' The summary.mpwgaim object has components that allow plotting of a
##' log-probability measures at each interval or marker including the
##' selected QTL which are labelled.  It is possible to obtain plots
##' of a subset the chromosomes or linkage groups.  By default, only
##' those chromosomes or linkage groups containing QTL (chr = 'qtl')
##' are plotted.  The completed plot can be saved as a postscript
##' file.
##' @title Plot a 'summary.mpwgaim' object
##' @param x An object of class 'summary.mpwgaim'
##' @param intervalObj An object of class 'interval'
##' @param chr A vector consisting of which chromosomes or linkage
##' groups to plot.  The default setting is 'qtl' which results in
##' plotting the chromosomes or linkage groups containing QTL.
##' @param plot Logical indicating whether or not to plot.  The
##' default is TRUE.
##' @param save.plot Logical indicating whether or not to save the
##' plot to a file.  Default is FALSE.
##' @param file File name of the saved plot file if save.plot is
##' TRUE.  If no name is given the file is called 'Rplot.ps' otherwise
##' it is the same supplied in this argument.
##' @param \ldots  Additional optional arguments passed to or from
##' other methods
##' @return  The components of the plot are saved in a list with
##' components 'df', a data frame with the chromosome, position and
##' log-probability used in the plot and the plot itself save as the
##' component 'Plot'.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##' @examples
##' \dontrun{
##' ## Simulated data example
##' require(asreml)
##' require(mpMap)
##' require(qtl)
##' map <- sim.map(len=rep(200,7), n.mar=rep(51,7), eq.spacing=TRUE, include.x=FALSE)
##' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
##' ## QTL to set up mpcross structure
##' qtl.mat <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##'                   2, 162, 0.354, -0.354, -0.354, -0.354,
##'                   5, 78, 0.354, -0.354, -0.354, 0.354),
##'                   nrow=3, ncol=6, byrow=TRUE)
##' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped,
##'                        qtl=qtl.mat, seed=5)
##' mpSim <- maporder(sim.dat)
##' mpInterval <- mpcross2int(mpSim, gen.type="mpInterval")
##' mpMarker <- mpcross2int(mpSim, gen.type="mpMarker")
##' ##  Model matrix for computing the contribution of the QTL
##' nqtl.col <- dim(sim.dat$qtlgeno$finals)[2]
##' mmat <- matrix(nrow=500, ncol=4*nqtl.col/2)
##' for (ii in 1:(nqtl.col/2)) {
##'   qtl.fac <- factor(sim.dat$qtlgeno$finals[,ii])
##'   mmat[,(4*ii-3):(4*ii)] <- model.matrix(~qtl.fac - 1)
##' }
##' ## Effects
##' qtl.sizes <- as.vector(t(qtl.mat[,3:6]))
##' qtl.effect <- mmat %*% qtl.sizes
##' ## Polygenic variance
##' pvar <- 0.5
##' ## Function to calculate approximate percentage variance for each QTL
##' perc.var <- function(qtl.mat, poly.var) {
##'     nfounders <- dim(qtl.mat)[2]-2
##'     prob <- 1/nfounders
##'     varq <- diag(rep(prob,nfounders)) - rep(prob,nfounders) %*% t(rep(prob,nfounders))
##'     gvar <- apply(qtl.mat[, -c(1,2)], 1, function(el, varq) sum((el %*% varq) * el), varq)
##'     totvar <- sum(gvar)+pvar
##'     perc.var <- 100*gvar/totvar
##'     round(perc.var,1)
##' }
##' ## Percentage variance for each QTL
##' percvar <- perc.var(qtl.mat, pvar)
##' percvar
##' ## Setup simulated data for analysis
##' ngeno <- 500
##' nrep <- 2
##' id <- factor(rep(paste0("L", 1:ngeno),nrep))
##' pheno.data <- data.frame(y = 0, Rep=nrep, id=id)
##' set.seed(10)
##' ee <- rnorm(ngeno*nrep,0,1)
##' uu <- rnorm(ngeno,0,1)
##' pheno.data$y <- 10 + rep(qtl.effect,2) + sqrt(1/2)*c(uu,uu) + ee
##' sim.asr0 <- asreml(y ~ 1, random = ~ id, data=pheno.data)
##' sim.qtl <- mpwgaim(sim.asr0, pheno.data, mpInterval, merge.by = "id",
##'                   verboseLev=0, gen.type="interval", na.method.X='include',
##'                   data.name = "sim.data")
##' sim.summ <- summary(sim.qtl, mpInterval)
##' sim.summ
##' plot(sim.summ, mpInterval)
##' ## Can repeat the analysis using mpMarker instead of mpInterval
##' }
##'
plot.summary.mpwgaim <- function(x, intervalObj, chr = "qtl",
                                 plot = TRUE, save.plot=FALSE, file = NULL, ...)
{
    dots <- list(...)
    cex <- dots$cex
    if (is.null(qtl <- x$qtl))
        stop("There are no QTLs to plot.")
    qtl <- gsub("X.", "", qtl)
    if (x$type == "interval")
        dist <- lapply(intervalObj$geno, function(el) {
            if (length(el$map) == 1) {
                el$dist <- 0.05
                names(el$dist) <- names(el$map)
                el$dist
            }
            else el$dist
        })
    else dist <- lapply(intervalObj$geno, function(el) {
        tel <- 0.05
        if (length(el$map) > 1)
            tel <- el$dist
        names(tel)[1] <- names(el$map)[1]
        tel
    })
    if (x$type == "interval")
        dist <- lapply(dist, function(el) (el[-1] + el[-length(el)])/2)
    ln <- length(qtl)
    y.lab <- "Log Probability: LOGP"
    oint <- unlist(x$LOGPtrace)
    lchn <- unlist(lapply(x$LOGPtrace, function(el) length(el)))
    chn <- rep(names(x$LOGPtrace), lchn)
    dint <- cbind.data.frame(LOGP = oint, chn = chn, dist = unlist(dist))
    which.qtl <- grep("X.", row.names(dint))
    dint$qtl <- rep("", dim(dint)[1])
    dint$qtl[which.qtl] <- "QTL"
    if (!is.null(chr)) {
        if(chr[1] == "qtl") {
            tmp <- strsplit(qtl, split="\\.")
            chr <- unique(unlist(lapply(tmp, function(el) el[1])))
        }
        dint <- dint[as.character(dint$chn) %in% chr, ]
        dint$chn <- factor(as.character(dint$chn))
    }
    prow <- length(unique(dint$chn))
    ppage <- 1
    if (prow > 5) {
        ppage <- prow %% 5 + 1
        prow <- 5
    }
    if(ppage > 1) devAskNewPage(ask=TRUE)
    if(save.plot) {
        if(is.null(file)) file <- "Rplot.ps"
        trellis.device(postscript, file=file, po=14)
    }
    ##
    myPlot <- xyplot(LOGP ~ dist | chn, type = "l", data = dint,
                     panel = function(x,y,subscripts,qtl) {
                         panel.xyplot(x,y, type="l")
                         panel.text(x,y+0.1,label=qtl[subscripts], cex=0.5)
                     }, qtl = dint$qtl,
                     ylab = y.lab, xlab = "Distance (cM)",
                     layout = c(1, prow), ylim = c(0, max(dint$LOGP)+0.5), ...)
    ##
    if(plot) print(myPlot)
    if(save.plot) {
        cat("Plots saved to ", file, "\n")
        dev.off()
    }
    if(ppage > 1) devAskNewPage(ask=FALSE)
    invisible(list(df = dint, Plot = myPlot))
}
