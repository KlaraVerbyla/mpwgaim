##  Method summary.mpwgaim
##
##' Prints a summary of the "'mpwgaim'" object in a presentable format
##'
##' After a QTL analysis using 'mpwgaim', a summary of the QTL can be
##' obtained using this function.  In addition the building blocks are
##' calculated for a plot of log-probabilities across the genome
##' including QTL.
##' @title Summary method for the class 'mpwgaim'
##' @param object an object of class 'mpwgaim'
##' @param intervalObj an object of class 'interval'
##' @param merge.by The genetic line variable used in merging
##' phenotypic and genetic data in the 'mpwgaim' fit.
##' @param \ldots Optional additional objects passed to or from other
##' methods
##' @return A object of class 'summary.mpwgaim, consisting of a list
##' with a data frame of QTL information to be used with the print
##' method, components 'LODtrace', 'qtl' and 'type' that are used for
##' the plotting method 'plot.summary.mpwgaim'
##' @seealso 'print.summary.mpwgaim', 'plot.summary.mpwgaim'
##' @author Ari Verbyla (ari.verbyla at csiro.au) and Klara Verbyla
##' (klara.verbyla at csiro.au)
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
summary.mpwgaim <- function(object, intervalObj, merge.by="id", ...)
{
##    require(asreml)  delete - ABZ
##    require(wgaim)  delete - ABZ
##    require(MASS)  delete - ABZ
    phenoData <- eval(object$call$data)
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (is.null(qtls <- object$QTL$qtl)) {
        cat("There are no significant putative QTLs\n")
        return()
    }
    else {
        sigma2 <- object$sigma2
        if (names(object$gammas.con[length(object$gammas.con)]) ==
            "Fixed")
            sigma2 <- 1
        nfounders <- object$QTL$nfounders
        LOGP <- LogProb <- c()
        Pvalue <- c()
        effects <- i.LOGP <- i.LogProb <- list()
        i.Pvalue <- list()
        for (ii in 1:length(qtls)) {
            qtl <- object$QTL$qtl[ii]
            whq <- unlist(eval(object$call$group[[qtl]]))
            n.int <- length(whq)
            gqtl <- gsub("X.", "", qtl) 
            cat(" Predict step for QTL ", gqtl, "\n")
            cat("===============================\n")
            qlev <- diag(n.int)
            qlist <- list(as.vector(qlev))
            names(qlist) <- qtl
            pred <- predict(object, classify = paste0('grp(',qtl,')'), only = paste0('grp(',qtl,')'),
                            levels = qlist, vcov = TRUE,
                            maxiter = 1)$predictions
            effects[[ii]] <- pred$pvals$predicted.value
            vcov <- pred$vcov
            d2.ii <- sum(effects[[ii]] %*% ginv(vcov) * effects[[ii]])
            i.d2.ii <- (effects[[ii]]^2) / diag(vcov)
            Pvalue[ii] <- pchisq(d2.ii, df = nfounders-1, lower.tail=FALSE)
            LOGP[ii] <- round(0.5 * log(exp(d2.ii), base = 10), 2)
            LogProb[ii] <- round(-log(Pvalue[ii], base=10), 2)
            i.Pvalue[[ii]] <- pchisq(i.d2.ii, df = 1, lower.tail=FALSE)/2
            i.LOGP[[ii]] <- round(0.5 * log(exp(i.d2.ii), base = 10), 3)
            i.LogProb[[ii]] <- round(-log(i.Pvalue[[ii]], base=10), 2)
        }
        iPvalues<-unlist(i.Pvalue)
        iLOGPs<-unlist(i.LOGP)
        iLogProbs <- unlist(i.LogProb)
        qtlmat <- mpgetQTL(object, intervalObj)
        Pvalues <- addons <- LOGPs <- LogProbs <- rep("", length(object$QTL$effects))
        ##
        seqs<-seq(from=1, to=length(Pvalues), by=nfounders)
        ##
        Pvalues[seqs] <- round(Pvalue, 3)
        LOGPs[seqs] <- round(LOGP, 2)
        LogProbs[seqs] <- round(LogProb, 2)
        ##
        qtlmat.expand <- matrix("", nrow=dim(qtlmat)[1]*nfounders, ncol=dim(qtlmat)[2])
        qtlmat.expand[seqs,] <- qtlmat
        wchr <- rep(qtlmat[,1], each=nfounders)
        wint <- rep(as.numeric(qtlmat[,2]), each=nfounders)
        ##
        gen.type <- intervalObj$gen.type
        if(gen.type == "mpPost") gen.type = "marker"
        if (gen.type == "mpInterval")
            gdat <- lapply(intervalObj$geno, function(el) el$intval)
        else
            gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
        if(intervalObj$gen.type != "mpPost") {
            if(any(lapply(gdat, function(x) dim(x)[2] %% nfounders) != 0))
                stop("Additive Genetic data is not consistent with the number of founders")
            nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
            lint <- unlist(lapply(nint, length))
            gnams <- paste("Chr", rep(names(intervalObj$geno), times = lint), unlist(nint), sep = ".")
            eindex <- rep(1:nfounders,length(gnams))
            gnams <- paste(rep(gnams,each=nfounders), eindex, sep=".")
        }
        if(intervalObj$gen.type == "mpPost") {
            gnams <- unlist(lapply(intervalObj$geno, function(el) dimnames(el$imputed.data)[[2]]))
        }
        geneticData <- do.call("cbind.data.frame", gdat)
        names(geneticData) <- gnams
        gD <- geneticData[, !(gnams %in% unlist(object$QTL$qtl.list))]
        gnam <- unlist(lapply(strsplit(names(object$gammas), "!"),
                              function(el) el[1]))
        gnam <- names(intervalObj$pheno)[names(intervalObj$pheno) %in%
                                         gnam]
        qtl.var <- c()
        for (ii in 1:length(object$QTL$qtl.list)) {
            probs <- geneticData[, object$QTL$qtl.list[[ii]]]
            av.prob <- apply(probs, 2, mean)
            prob.mat <- av.prob %*% t(av.prob)
            qtl.var[ii] <- sum((effects[[ii]] %*% (diag(av.prob) - prob.mat)) * effects[[ii]])
        }
        coeff.other <- sum(apply(gD, 2, mean)^2)
        other.var <- sigma2 * object$gammas[grep(paste("!", gnam,
                                                       sep = ""), names(object$gammas))]
######  next command is needed if there is no replication of lines and hence no id variance
        if(length(other.var) == 0) other.var <- sigma2
        coeff.other <- c(coeff.other, rep(1, length(other.var)))
        other.var <- c(sigma2 * object$gammas[grep("ints",
                                                   names(object$gammas))]/10000, other.var)
        nonqtl.var <- coeff.other*other.var
        full.var <- sum(qtl.var) + sum(nonqtl.var)
        perc.var <- 100 * qtl.var/full.var
        names(perc.var) <- qtls
        adds <- round(perc.var, 1)
        addons[seqs] <- adds
        effects <- unlist(effects)
        qtlmat <- cbind(qtlmat.expand, rep(1:nfounders, length(qtls)),
                                   round(effects, 3), round(iPvalues, 3),
                        round(iLogProbs,2), Pvalues, addons, LogProbs)
        if (object$QTL$type == "interval") {
            collab <- c("Chromosome", "Interval No.", "Left Marker", "dist (cM)",
                        "Right Marker", "dist (cM)", "Founder", "Size",
                        "Founder Prob", "Founder LOGP", "Prob", "% var", "LOGP")
            fndrpos <- 7
            qtlmat <- qtlmat[order(wchr, wint, as.numeric(qtlmat[, fndrpos])), ]
        }
        else {
            collab <- c("Chromosome", "Marker No.", "Marker", "dist (cM)",
                        "Founder", "Size", "Founder Prob", "Founder LOGP",
                        "Prob", "% var", "LOGP")
            fndrpos <- 5
            qtlmat <- qtlmat[order(wchr, wint, as.numeric(qtlmat[, fndrpos])), ]
        }
        dimnames(qtlmat) <- list(as.character(1:length(effects)), collab)
        Founder <- rep(intervalObj$founders, length(qtls))
        qtlmat[, "Founder"] <- Founder
        qtlmat <- qtlmat[,-2]
        qtlmat <- as.data.frame(qtlmat)
        rownames(qtlmat) <- 1:dim(qtlmat)[1]
    }
    ##
################ Components for plot
    ##
    whq <- unlist(eval(object$call$group$ints))
    n.int <- length(whq)
    cat(" Predict step for non-QTL loci \n")
    cat("===============================\n")
    qlev <- diag(n.int)
    qlist <- list(as.vector(qlev))
    names(qlist) <- "ints"
    pv <- predict(object, classify = "grp(\"ints\")", only = "grp(\"ints\")",
                  levels = qlist, vcov = TRUE, maxiter = 1)
    gamma <- pv$gammas[grep("ints", names(object$gammas))]
    atilde <- pv$predictions$pvals[, "predicted.value"]
    pev <- pv$predictions$vcov
    vatilde <- sigma2 * gamma * diag(n.int) - pev
    p <- dim(gD)[2]
    other <- intervalObj$pheno[, merge.by]
    whg <- !duplicated(phenoData[, merge.by])
    whg <- other %in% phenoData[whg, merge.by]
    ids <- as.character(other[whg])
    q <- length(ids)
    gDnams <- dimnames(gD)[[2]]
    if (p > q) {
        mats <- gD[whg, ]
        tmat <- t(mats)
        xsvd <- svd(crossprod(tmat))
        xsvd.half <- xsvd.half <- t(xsvd$v %*% (t(xsvd$u) * sqrt(xsvd$d)))
        xsvd.inv <- solve(xsvd.half)
        xtra <- tmat %*% xsvd.inv
        dimnames(xtra) <- list(gDnams, 1:dim(xtra)[2])
        atilde <- xtra %*% atilde
        names(atilde) <- gDnams
    }
    else {
        dimnames(vatilde) <- list(gDnams, gDnams)
    }
    ##
    split.gDnams <- strsplit(gDnams, split='\\.')
    mgnams <- unique(unlist(lapply(split.gDnams,
                                   function(x) paste(x[1], x[2], x[3], sep='.'))))
    LogProb.ii <- c()
    noint <- length(mgnams)
    for (ni in 1:noint) {
        wh.m <- gDnams[grep(paste(mgnams[ni], "\\.", sep=""), gDnams)]
        if(p > q) {
            vcov.ii <- xtra[wh.m,] %*% vatilde %*% t(xtra[wh.m,])
        }
        else {
            vcov.ii <- vatilde[wh.m, wh.m]
        }
        d2.ii <- t(atilde[wh.m]) %*% ginv(vcov.ii) %*% atilde[wh.m]
        Pvalue.ii <- pchisq(d2.ii, df = nfounders-1, lower.tail=FALSE)
        LogProb.ii[ni] <- round(-log(Pvalue.ii, base=10), 2)
    }
    names(LogProb.ii) <- mgnams
    LogProb.ii[is.na(LogProb.ii)] <- 0
    names(LogProb) <- qtls
    allLogProb <- vector(mode = "list", length=length(intervalObj$geno))
    names(allLogProb) <- names(intervalObj$geno)
    all.tmp <- c(LogProb, LogProb.ii)
    split.tmp <- strsplit(names(all.tmp), split='\\.')
    wchr <- unlist(lapply(split.tmp, function(el) el[2]))
    wint <- as.numeric(unlist(lapply(split.tmp, function(el) el[3])))
    for (ii in names(allLogProb)) {
        allLogProb[[ii]] <- all.tmp[wchr == ii]
        ord.ii <- order(wint[wchr == ii])
        allLogProb[[ii]] <- allLogProb[[ii]][ord.ii]
    }
    outObject <- list(summary=qtlmat, LOGPtrace=allLogProb, qtl = qtls,
                      type = object$QTL$type)
    class(outObject) <- "summary.mpwgaim"

    invisible(outObject)
}