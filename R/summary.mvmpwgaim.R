##  Method summary.mvmpwgaim
##
##' Provide a summary of the QTL found using mvmpwgaim.
##'
##' A summary data frame is created for multi-trait or
##' multi-environment QTL analysis using \code{mvmpwgaim}.
##' @title Summary method for the class 'mvmpwgaim'
##' @param object An object of class \code{mvmpwgaim} after fitting
##' using \code{mvmpwgaim.asreml}
##' @param intervalObj An interval object used in the fit of the
##' model resulting in object.
##' @param by The name of the genotype factor for the genetic lines
##' in the analysis
##' @param LOGP Logical set to \code{TRUE}, so that -log(Probability)
##' is included in the summary of the QTL found.
##' @param \ldots Optional additional arguments.
##' @return A data frame with a summary of the QTL found, including
##' location (chromosome or linkage group, the interval or marker
##' including genetic distance), size of the QTL for each founder by
##' trait or environment combination, a probability measure of the
##' 'significance" of the QTL, percentage variance for each QTL for
##' each trait or environment and a LOGP measure if indicated by the
##' \code{LOGP} argument.
##' @author Ari Verbyla (ari.verbyla at csiro.au) and Klara Verbyla
##' (klara.verbyla at csiro.au)
##' @export
##' @examples
##' \dontrun{
##' ## Simulation of genetic and phenotypic data for analysis
##' require(asreml)
##' require(mpMap)
##' require(qtl)
##' ## Generate a linkage map
##' map <- sim.map(len=rep(200,7), n.mar=rep(51,7), eq.spacing=TRUE, include.x=FALSE)
##' ## Set up a pedigree  for a 4-way cross with 500 RILs and 6 generations of selfing
##' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
##' ## Dummy QTL structure
##' qtl.mat <- matrix(data=c(1, 142, 0.354, -0.354, -0.354,  0.354,
##'                          2, 162, 0.354, -0.354, -0.354, -0.354,
##'                          5, 78, 0.354, -0.354, -0.354, 0.354),
##'                          nrow=3, ncol=6, byrow=TRUE)
##' ## Simulate mpcross object
##' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped,
##'                        qtl=qtl.mat, seed=5)
##' mpSim <- maporder(sim.dat)
##' ## Create interval objects (see mpcross2int)
##' mpInterval <- mpcross2int(mpSim, gen.type="mpInterval")
##' mpMarker <- mpcross2int(mpSim, gen.type="mpMarker")
##' ##  Model matrix for computing the contribution of the QTL
##' nqtl.col <- dim(sim.dat$qtlgeno$finals)[2]
##' mmat <- matrix(nrow=500, ncol=4*nqtl.col/2)
##' for (ii in 1:(nqtl.col/2)) {
##' qtl.fac <- factor(sim.dat$qtlgeno$finals[,ii])
##' mmat[,(4*ii-3):(4*ii)] <- model.matrix(~qtl.fac - 1)
##' }
##' ##  Effects for each environment
##' ## Environment 1
##' qtl.mat1 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##'                           2, 162, 0.354, -0.354, -0.354, 0.354,
##'                           5, 78, 0.354, -0.354, -0.354, 0.354),
##'                           nrow=3, ncol=6, byrow=TRUE)
##' qtl.sizes1 <- as.vector(t(qtl.mat1[,3:6]))
##' qtl.effect1 <- mmat %*% qtl.sizes1
##' ## Environment 2
##' qtl.mat2 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##'                           2, 162, 0.354, -0.354, -0.354, 0.354,
##'                           5, 78, 0, 0, 0, 0),
##'                           nrow=3, ncol=6, byrow=TRUE)
##' qtl.sizes2 <- as.vector(t(qtl.mat2[,3:6]))
##' qtl.effect2 <- mmat %*% qtl.sizes2
##' ## Environment 3
##' qtl.mat3 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##'                           2, 162, -0.354, 0.354, 0.354, -0.354,
##'                           5, 78, 0.354, -0.354, -0.354, 0.354),
##'                           nrow=3, ncol=6, byrow=TRUE)
##' qtl.sizes3 <- as.vector(t(qtl.mat3[,3:6]))
##' qtl.effect3 <- mmat %*% qtl.sizes3
##' ## Polygenic variance
##' pvar_0.5
##' ## Function to calculate approximate percentage variance for each QTL
##' perc.var <- function(qtl.mat, poly.var) {
##'    nfounders <- dim(qtl.mat)[2]-2
##'    prob <- 1/nfounders
##'    varq <- diag(rep(prob,nfounders)) - rep(prob,nfounders) %*% t(rep(prob,nfounders))
##'    gvar <- apply(qtl.mat[, -c(1,2)], 1, function(el, varq) sum((el %*% varq) * el), varq)
##'    totvar <- sum(gvar)+pvar
##'    perc.var <- 100*gvar/totvar
##'    round(perc.var,1)
##' }
##' ## Percentage variance for each QTL in each environment
##' ## Environment 1
##' percvar1 <- perc.var(qtl.mat1, pvar)
##' percvar1
##' ## Environment 2
##' percvar2 <- perc.var(qtl.mat2, pvar)
##' percvar2
##' ## Environment 3
##' percvar3 <- perc.var(qtl.mat3, pvar)
##' percvar3
##' ## Setup simulated data for analysis
##' ntrait <- 3
##' ngeno <- 500
##' nrep <- 2
##' Trait <- factor(rep(rep(1:ntrait, ngeno), nrep))
##' gvar <- matrix(c(1.0,0.9,0.7,
##'                  0.9,1.0,0.7,
##'                  0.7,0.7,1.0), ncol=3)
##' gchol <- t(chol(gvar))
##' id <- factor(rep(rep(paste0("L", 1:ngeno),each=ntrait),nrep))
##' pheno.data <- data.frame(y = 0, Trait=Trait, Rep=nrep, id=id)
##' qtleffects <- as.matrix(0, ncol=1, nrow=ntrait*ngeno)
##' qtleffects[rep(1:3,ngeno)==1] <- qtl.effect1
##' qtleffects[rep(1:3,ngeno)==2] <- qtl.effect2
##' qtleffects[rep(1:3,ngeno)==3] <- qtl.effect3
##' set.seed(10)
##' ee <- rnorm(ntrait*ngeno*nrep,0,1)
##' geno <- rnorm(ntrait*ngeno,0,1)
##' uu <- as.vector(gchol%*%matrix(geno, nrow=ntrait))
##' pheno.data$y <- rep(c(9,10,12), ngeno*nrep)+
##'     qtleffects + sqrt(1/2)*c(uu,uu) + ee
##' ord <- order(pheno.data$Trait)
##' pheno.data<-pheno.data[ord,]
##' ## Multi-environment analysis
##' ## Multi-environment preliminary models
##' test.asr0 <- asreml(y ~ Trait - 1, random = ~ diag(Trait):id, data=pheno.data)
##' test.asr1 <- asreml(y ~ Trait - 1, random = ~ fa(Trait,1):id, data=pheno.data)
##' ##  MVMPWGAIM analysis
##' test.qtl <- mvmpwgaim(test.asr0, test.asr1, pheno.data, mpInterval, merge.by = "id",
##'                       Trait="Trait", verboseLev=0, main.effects = TRUE,
##'                       data.name = "mvy.data", gen.type="interval", na.method.X='include',
##'                       workspace=4e7, pworkspace=4e7)
##' ## Check for Trait by QTL interactions and fit a reduced model if needed
##' test.aqtl <- remlrt.mvmpwgaim(test.qtl, mvy.data, TypeI=0.05)
##' ##  Summary for MVMPWGAIM analysis
##' test.summ <- summary(test.aqtl, mpInterval)
##' test.summ
##' }
##'
summary.mvmpwgaim <-
function (object, intervalObj, by = "id", LOGP = TRUE, ...)
{
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (is.null(object$QTL$qtl)) {
        cat("There are no significant putative QTL's\n")
        return(character(0))
    }
    else {
        phenoData <- eval(object$call$data)
        sc <- names(object$gammas.con)[length(object$gammas.con)]
        if (sc == "Fixed")
            sigma2 <- 1
        else sigma2 <- object$sigma2
        trait <- object$QTL$trait
        n.trait <- length(object$QTL$trait.levels)
        nfounders <- object$QTL$nfounders
        qtl <- unlist(object$QTL$qtl)
        n.qtl <- length(qtl)
        if (object$QTL$main)
            qtl.type <- object$QTL$qtl.type
        else qtl.type <- rep("interaction", n.qtl)
        qtl.m <- qtl[qtl.type == "main"]
        qtl.i <- qtl[qtl.type != "main"]
        grp1 <- grep("X\\.", object$factor.names)
        which1 <- object$factor.names[grp1]
        grp2 <- grep(object$QTL$trait, which1)
        if (length(grp2) > 0) {
            which2 <- which1[grp2]
            which2.str <- strsplit(which2, split = ":")
            which3 <- c()
            for (ii in 1:length(which2.str)) {
                which.ii <- grep("X\\.", which2.str[[ii]])
                which3[ii] <- which2.str[[ii]][which.ii]
            }
        }
        else which3 <- c()
        effects <- Pvalue <- LogProb <- i.Pvalue <- i.LogProb <- Trait.list <- qtl.int <- founder.list <- list()
        if (length(qtl.m) > 0) {
            qtl.namm <- gsub("X.", "", qtl.m)
            for (ii in 1:length(qtl.m)) {
                qtl.ii <- qtl.m[ii]
                whq <- unlist(eval(object$call$group[[qtl.ii]]))
                n.int <- length(whq)
                cat(" Predict step for QTL ", qtl.namm[ii], " for summary \n")
                cat("==========================================\n")
                qlev <- diag(n.int)
                qlist <- list(as.vector(qlev))
                names(qlist) <- qtl.ii
                pred <- predict(object, classify = paste0("grp(",
                  qtl.ii, ")"), only = paste0("grp(", qtl.ii,
                  ")"), levels = qlist, vcov = TRUE, maxiter = 1,
                  data = phenoData)$predictions
                effects[[ii]] <- pred$pvals$predicted.value
                names(effects)[ii] <- qtl.ii
                vcov <- pred$vcov
                d2.ii <- sum(effects[[ii]] %*% ginv(vcov) * effects[[ii]])
                Pvalue[[ii]] <- rep(pchisq(d2.ii, df = nfounders -
                  1, lower.tail = FALSE), n.trait)
                LogProb[[ii]] <- round(-log(Pvalue[[ii]], base = 10),
                  2)
                i.d2.ii <- rep((effects[[ii]]^2)/diag(vcov),
                  n.trait)
                i.Pvalue[[ii]] <- pchisq(i.d2.ii, df = 1, lower.tail = FALSE)/2
                i.LogProb[[ii]] <- round(-log(i.Pvalue[[ii]],
                  base = 10), 2)
                Trait.list[[ii]] <- rep(object$QTL$trait.levels,
                  each = nfounders)
                founder.list[[ii]] <- rep(pred$pvals[, qtl.ii],
                  n.trait)
                effects[[ii]] <- rep(effects[[ii]], n.trait)
                tmp.split <- strsplit(qtl.ii, split = "\\.")
                wchr <- unlist(lapply(tmp.split, function(el) el[2]))
                wint <- as.numeric(unlist(lapply(tmp.split, function(el) el[3])))
                qtl.int[[ii]] <- rep(paste(wchr, wint, sep = "."),
                  nfounders * n.trait)
            }
        }
        n.main <- length(qtl.m)
        if (length(which3) > 0) {
            qtl.nami <- gsub("X.", "", qtl.i)
            for (ii in 1:length(which3)) {
                jj <- ii + n.main
                whq <- unlist(eval(object$call$group[qtl.i[ii]]))
                n.int <- length(whq)
                cat(" Predict step for QTL ", qtl.nami[ii], " for summary \n")
                cat("==========================================\n")
                qlev <- diag(n.int)
                qlist <- list(as.vector(qlev))
                names(qlist) <- znt.me.ii <- which3[ii]
                if (!object$QTL$main)
                  znt.me.ii <- NULL
                pred.ii <- predict(object, classify = which2[ii],
                  only = c(which2[ii], znt.me.ii), levels = qlist,
                  maxiter = 1, vcov = TRUE, data = phenoData)$predictions
                effects[[jj]] <- pred.ii$pvals$predicted.value
                names(effects)[jj] <- qtl.i[ii]
                vcov <- pred.ii$vcov
                Pvalue[[jj]] <- 0
                for (kk in 1:n.trait) {
                  kk.left <- (kk - 1) * nfounders + 1
                  kk.right <- kk * nfounders
                  d2.ii <- sum(effects[[jj]][kk.left:kk.right] %*%
                    ginv(vcov[kk.left:kk.right, kk.left:kk.right]) *
                    effects[[jj]][kk.left:kk.right])
                  Pvalue[[jj]] <- c(Pvalue[[jj]], pchisq(d2.ii,
                    df = nfounders - 1, lower.tail = FALSE))
                }
                Pvalue[[jj]] <- Pvalue[[jj]][-1]
                LogProb[[jj]] <- round(-log(Pvalue[[jj]], base = 10),
                  2)
                i.d2.ii <- (effects[[jj]]^2)/diag(vcov)
                i.Pvalue[[jj]] <- pchisq(i.d2.ii, df = 1, lower.tail = FALSE)/2
                i.LogProb[[jj]] <- round(-log(i.Pvalue[[jj]],
                  base = 10), 2)
                Trait.list[[jj]] <- object$QTL$trait.levels[as.numeric(pred.ii$pvals[,
                  trait])]
                founder.list[[jj]] <- pred.ii$pvals[, qtl.i[ii]]
                tmp.split <- strsplit(qtl.i[ii], split = "\\.")
                wchr <- unlist(lapply(tmp.split, function(el) el[2]))
                wint <- as.numeric(unlist(lapply(tmp.split, function(el) el[3])))
                qtl.int[[jj]] <- rep(paste(wchr, wint, sep = "."),
                  nfounders * n.trait)
            }
        }
        tmp.split <- strsplit(unlist(qtl.int), split = "\\.")
        wchr <- unlist(lapply(tmp.split, function(el) el[1]))
        wint <- as.numeric(unlist(lapply(tmp.split, function(el) el[2])))
        qtlmat <- data.frame(int = unlist(qtl.int), wchr = wchr,
            wint = wint, unlist(Trait.list), Founder = unlist(founder.list),
            Size = unlist(effects), Founder.Prob = unlist(i.Pvalue),
            Founder.LOGP = unlist(i.LogProb))
        names(qtlmat)[4] <- trait
        seqs <- seq(from = 1, to = dim(qtlmat)[1], by = nfounders)
        qtlmat$Prob <- NA
        qtlmat$Prob[seqs] <- unlist(Pvalue)
        qtlmat$LOGP <- NA
        qtlmat$LOGP[seqs] <- unlist(LogProb)
        gen.type <- intervalObj$gen.type
        if (gen.type == "mpInterval")
            gdat <- lapply(intervalObj$geno, function(el) el$intval)
        else gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
        if (any(lapply(gdat, function(x) dim(x)[2]%%nfounders) !=
            0))
            stop("Additive Genetic data is not consistent with the number of founders")
        nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
        lint <- unlist(lapply(nint, length))
        gnams <- paste("Chr", rep(names(intervalObj$geno), times = lint),
            unlist(nint), sep = ".")
        eindex <- rep(1:nfounders, length(gnams))
        gnams <- paste(rep(gnams, each = nfounders), eindex,
            sep = ".")
        geneticData <- do.call("cbind.data.frame", gdat)
        names(geneticData) <- gnams
        gD <- geneticData[, !(gnams %in% unlist(object$QTL$qtl.list))]
        qtl.var <- matrix(nrow = n.trait, ncol = length(object$QTL$qtl.list))
        for (ii in 1:length(object$QTL$qtl.list)) {
            probs <- geneticData[, object$QTL$qtl.list[[ii]]]
            av.prob <- apply(probs, 2, mean)
            prob.mat <- av.prob %*% t(av.prob)
            for (kk in 1:n.trait) {
                kk.left <- (kk - 1) * nfounders + 1
                kk.right <- kk * nfounders
                qtl.var[kk, ii] <- sum((effects[[ii]][kk.left:kk.right] %*%
                  (diag(av.prob) - prob.mat)) * effects[[ii]][kk.left:kk.right])
            }
        }
        coeff.other <- sum(apply(gD, 2, mean)^2)
        poly <- sigma2 * object$gammas[grep(paste(by, "!", sep = ""),
            names(object$gammas))]
        ints <- sigma2 * object$gammas[grep("ints", names(object$gammas))]/10000
        if (n.trait > 2) {
            n.fa <- length(poly)/n.trait - 1
            p.psi <- poly[1:n.trait]
            p.lambda <- matrix(poly[(n.trait + 1):length(poly)],
                ncol = n.fa, byrow = FALSE)
            p.var <- diag(p.lambda %*% t(p.lambda) + diag(p.psi))
            n.fa <- length(ints)/n.trait - 1
            i.psi <- ints[1:n.trait]
            i.lambda <- matrix(ints[(n.trait + 1):length(ints)],
                ncol = n.fa, byrow = FALSE)
            i.var <- diag(i.lambda %*% t(i.lambda) + diag(i.psi))
        }
        if (n.trait == 2) {
            p.var <- i.var <- matrix(0, 2, 2)
            p.var[1, 1] <- poly[2]
            p.var[2, 2] <- poly[3]
            p.var[1, 2] <- p.var[2, 1] <- poly[1] * sqrt(poly[2] *
                poly[3])
            p.var <- diag(p.var)
            i.var[1, 1] <- ints[2]
            i.var[2, 2] <- ints[3]
            i.var[1, 2] <- i.var[2, 1] <- ints[1] * sqrt(ints[2] *
                ints[3])
            i.var <- diag(i.var)
        }
        allVars <- cbind(qtl.var, i.var, p.var)
        perc.var <- as.vector(allVars[, 1:n.qtl] * 100/rowSums(allVars))
        qtlmat$"% var" <- NA
        qtlmat$"% var"[seqs] <- perc.var
        qtlmat[, "Size"] <- round(qtlmat[, "Size"], 3)
        qtlmat[, "Prob"] <- round(qtlmat[, "Prob"], 4)
        qtlmat[, "Founder.Prob"] <- round(qtlmat[, "Founder.Prob"],
            4)
        qtlmat[, "% var"] <- round(qtlmat[, "% var"], 1)
        qtlmat[, "LOGP"] <- round(qtlmat[, "LOGP"], 1)
        qtlmat <- qtlmat[order(qtlmat[, "wchr"], qtlmat[, "wint"],
            qtlmat[, trait], qtlmat[, "Founder"]), ]
        qtlmat$Order <- 1:dim(qtlmat)[1]
        qtlmat$wchr <- qtlmat$wint <- NULL
        tmp.split <- strsplit(qtl, split = "\\.")
        wchr <- unlist(lapply(tmp.split, function(el) el[2]))
        wint <- as.numeric(unlist(lapply(tmp.split, function(el) el[3])))
        intmat <- mpwgaim:::getmvQTL(object, intervalObj)
        intmat <- cbind(paste(wchr, wint, sep = "."), intmat)
        intmat <- as.data.frame(intmat)[, -3]
        names(intmat)[1] <- "int"
        qtl.table <- merge(intmat, qtlmat, by = "int")
        qtl.table <- qtl.table[, -1]
        if (object$QTL$type == "interval") {
            names(qtl.table)[1:5] <- c("Chromosome", "Left Marker",
                "dist (cM)", "Right Marker", "dist (cM)")
            traitpos <- 6
            fndrpos <- 7
        }
        else {
            names(qtl.table)[1:3] <- c("Chromosome", "Marker",
                "dist (cM)")
            traitpos <- 4
            fndrpos <- 5
        }
        qtl.table <- qtl.table[order(qtl.table[, "Order"]), ]
        qtl.table[, "Order"] <- NULL
        qtl.table[, fndrpos] <- intervalObj$founders[qtl.table[,
                                                               fndrpos]]
        oseqs <- c(1:(n.trait * nfounders * n.qtl))[!(1:(n.trait *
            nfounders * n.qtl) %in% seq(from = 1, to = (n.trait *
            nfounders * n.qtl), by = nfounders))]
        oseqs2 <- c(1:(n.trait * nfounders * n.qtl))[!(1:(n.trait *
            nfounders * n.qtl) %in% seq(from = 1, to = (n.trait *
            nfounders * n.qtl), by = (nfounders * n.trait)))]
        final.table <- as.matrix(qtl.table)
        final.table[is.na(final.table)] <- ""
        final.table[oseqs, traitpos] <- ""
        final.table[oseqs2, 1:(traitpos - 1)] <- ""
        final.table <- as.data.frame(final.table)
        rownames(final.table) <- 1:dim(final.table)[1]
    }
    invisible(final.table)
}
