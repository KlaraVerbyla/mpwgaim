## Method mpwgaim.asreml
##
##' Fits an iterative Multi-Parent Whole Genome Average Interval
##' Mapping (mpwgaim) model for QTL detection
##'
##' This function implements the whole genome average interval mapping
##' approach for multi-parent populations.  It includes code to handle
##' high-dimensional problems.  The size of computations in model
##' fitting is the number of genetic lines rather than the number of
##' markers or intervals.  The use of mpwgaim requires a base model to
##' be fitted which includes all genetic and non-genetic terms that
##' are relevant for an analysis without marker information.  The
##' marker information must be contained in an interval object
##' ('intervalObj') and these effects are added to the phenotypic data
##' ('phenoData') to allow fitting of models with QTL.  Note that
##' analysis using 'markers' or 'intervals' can be carried out
##' ('gen.type' is set by default to 'interval').  All QTL selected
##' are random effects in the final model.
##' @title mpwgaim method for class \code{"asreml"}
##' @param baseModel a model object of class "'asreml'" usually
##' representing a base model with which to build the qtl model.
##' @param phenoData a data frame containing the phenotypic elements
##' used to fit 'baseModel'. This data is checked against the data
##' used in fitting the base model
##' @param intervalObj a list object containing the genotypic data,
##' usually an "'interval'" object obtained from using
##' 'mpcross2int'. This object may contain many more markers than
##' observations (see Details).
##' @param merge.by a character string or name of the column(s) in
##' 'phenoData' and 'intervalObj' to merge the phenotypic and
##' genotypic data sets.
##' @param gen.type a character string determining the type of
##' genetic data to be used in the analysis. Possibilities are
##' "'interval'" and "'markers'". The default is "'interval'". (see
##' Details).
##' @param TypeI a numerical value determining the level of
##' significance for detecting a QTL. The default is 0.05.
##' @param attempts An integer representing the number of attempts at
##' convergence for the QTL model. The default is 5.
##' @param data.name character string that represents the name of the
##' data frame for the final model fit using mpwgaim.  If no data.name
##' is specified, the file is saved in the working directory as the
##' name of the response dot data.
##' @param trace An automatic tracing facility. If 'trace = TRUE'
##' then all 'asreml' output is piped to the screen during the
##' analysis. If 'trace = "file.txt"', then output from all asreml
##' models is piped to "'file.txt'". Both trace mechanisms will
##' display a message if a QTL is detected.
##' @param verboseLev numerical value, either 0 or 1, determining the
##' level of tracing outputted during execution of the algorithm.  A 0
##' value will produce the standard model fitting output from the
##' fitted ASReml models involved in the forward selection. A value of
##' 1 will add a table of interval outlier statistics for each
##' iteration.
##' @param \ldots Any other extra arguments to be passed to each of
##' the 'asreml' calls. These may also include 'asreml.control'
##' arguments.
##' @return An object of class "'mpwgaim'" which also inherits the
##' class "'asreml'" by default. The object returned is actually an
##' 'asreml' object (see 'asreml.object') with the addition of
##' components from the QTL detection listed below.  \item{QTL:}{ A
##' list of components from the significant QTL detected including a
##' character vector of the significant QTL along with a vector of the
##' QTL effect sizes. There are also a number of diagnostics that can
##' be found in the 'diag' component.}
##' @seealso 'print.mpwgaim', 'summary.mpwgaim, 'plot.summary.mpwgaim'
##' @author Ari Verbyla (ari.verbyla at csiro.au) and Klara Verbyla
##' (klara.verbyla at csiro.au)
##' @references Verbyla, A. P., George, A. W., Cavanagh, C. C. and
##' Verbyla, K. L. (2014).  Whole genome QTL analysis for MAGIC.
##' Theoretical and Applied Genetics.  DOI:10.1007/s00122-014-2337-4.
##' @export
##' @examples ## Simulated data example
##' \dontrun{
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
##' ## Use Hidden Markov model for probability calculations
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
mpwgaim.asreml <- function(baseModel, phenoData, intervalObj, merge.by = NULL,
                           gen.type = "interval", TypeI = 0.05, attempts = 5,
                           data.name = NULL, trace = TRUE, verboseLev = 0, ...)
{
    if (missing(phenoData))
        stop("phenoData is a required argument.")
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (inherits(intervalObj, "mpInterval"))
        nfounders <- intervalObj$nfounders
    if (is.null(nfounders))
        stop("intervalObj doe not have founder information")
##    require(asreml) delete - ABZ
##    require(wgaim) delete - ABZ
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(other <- intervalObj$pheno[[merge.by]]))
        stop("Genotypic data does not contain column \"", merge.by,
            "\".")
    if (is.null(phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
            "\".")
    mby <- pmatch(as.character(other), as.character(phenoData[[merge.by]]))
    if (all(is.na(mby)))
        stop("Names in Genotypic \"", merge.by, "\" column do not match any
names in Phenotypic \"",
            merge.by, "\" column.")
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    if (gen.type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$intval)
    else gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    if(any(lapply(gdat, function(x) dim(x)[2] %% nfounders) != 0))
        stop(" Genetic data is not consistent with the number of founders")
    nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
    lint <- unlist(lapply(nint, length))
    gnams <- paste("Chr", rep(names(intervalObj$geno), times = lint),
                   unlist(nint), sep = ".")
    eindex <- rep(1:nfounders,length(gnams))
    gnams <- paste(rep(gnams,each=nfounders), eindex, sep=".")
    geneticData <- cbind.data.frame(other, do.call("cbind", gdat))
    names(geneticData) <- c(merge.by, gnams)
    dnams <- names(phenoData)
    basedata <- eval(baseModel$call$data)
    bnams <- names(basedata)
    if (any(is.na(pmatch(bnams, dnams))))
        stop("Some baseModel data names do not match phenoData names")
    whn <- unlist(lapply(basedata, is.numeric))
    whn <- whn[whn][1]
    diff <- unique(abs(basedata[[names(whn)]] - phenoData[[names(whn)]]))
    if (length(diff[!is.na(diff)]) > 1)
        stop("Phenotypic data is not in same order as baseModel data.\n Try reordering
phenoData apppropriately\n")
    int.cnt <- 2:dim(geneticData)[2]
    mD <- mpmergeData(phenoData, geneticData, merge.by)
    cnt <- mD$cnt
    asdata <- mD$asdata
    state <- rep(1, length(int.cnt))
    names(state) <- names(geneticData[, int.cnt])
    baseModel$call$data <- quote(asdata)
    baseModel$QTL$qtl
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
    }
    add.qtl <- baseModel
    add.qtl$call$group$ints <- cnt
    if(!is.null(baseModel$call$random)){add.form <- as.formula(paste("~ idv(grp('ints')) + ."))}
    if(is.null(baseModel$call$random)){add.form <- as.formula(paste("~ idv(grp('ints'))"))}
    class(add.qtl) <- "wgaim"
    cat("\nRandom Effects QTL Model Iteration (1):\n")
    cat("=========================================\n")
    add.qtl$call$data <- quote(asdata)
    add.qtl <- updateMPWgaim(add.qtl, asdata, attempts = attempts, random. = add.form, ...)
    update <- FALSE
    which.i <- 1
    ##
    dmat <- data.frame(L0 = 0, L1 = 0, Statistic = 0, Pvalue = 0)
    vl <- cl <- oint <- blups <- list()
    qtl.list <- list()
    qtl <- coeff.var <- c()
    ##
    repeat {
        if (update) {
            add.qtl$call$group[[last.qtl]] <- baseModel$call$group[[last.qtl]] <- qtls.cnt
            add.qtl$call$group$ints <- cnt
            if(!is.null(baseModel$call$random)){ qtl.form <- as.formula(paste("~ . + idv(grp(",last.qtl,"))", sep=""))}
            if(is.null(baseModel$call$random)){ qtl.form <- as.formula(paste("~ idv(grp(",last.qtl,"))", sep=""))}
            cat("\nRandom Effects QTL Model Iteration (", which.i,
                "):\n")
            cat("========================================\n")
            baseModel$call$data <- quote(asdata)
            baseModel <- mp.v.modify(baseModel)
            baseModel <- updateMPWgaim(baseModel, asdata, attempts = attempts,
                                     random. = qtl.form, ...)
            cat("\nRandom Effects QTL plus Interval/Marker Model Iteration (",
                which.i, "):\n")
            cat("=============================================================\n")
            add.qtl$call$data <- quote(asdata)
            if(!is.null(add.qtl$call$random)){ qtl.form <- as.formula(paste("~ . + idv(grp(",last.qtl,"))", sep=""))}
            if(is.null(add.qtl$call$random)){ qtl.form <- as.formula(paste("~ idv(grp(",last.qtl,"))", sep=""))}
            add.qtl <- mp.v.modify(add.qtl)
            add.qtl <- updateMPWgaim(add.qtl, asdata, attempts = attempts,
                                   random. = qtl.form, ...)
            list.coefs <- add.qtl$coefficients$random
            zind <- grep("X\\.", names(list.coefs))
            list.coefs <- list.coefs[zind]
            cl[[which.i - 1]] <- list.coefs
            vl[[which.i - 1]] <- add.qtl$vcoeff$random[zind]
        }
        baseLogL <- baseModel$loglik
        stat <- 2 * (add.qtl$loglik - baseLogL)
        pvalue <- (1 - pchisq(stat, 1))/2
        if(stat < 0) {
            stat <- 0
            pvalue <- 1
        }
        if(pvalue < 0) pvalue <- 0
        cat("\nLikelihood Ratio Test Statistic: ", stat, ", P-value: ", pvalue,"\n")
        dmat[which.i, ] <- c(baseLogL, add.qtl$loglik, stat, pvalue)
        if (pvalue > TypeI)
            break
        if (mD$p > mD$q)
            add.qtl$xtra <- mD$xtra
        add.qtl$call$data <- quote(asdata)
        pick <- mpqtl.pick(add.qtl, intervalObj, asdata, gen.type, state, verboseLev)
        state <- pick$state
        qtl.list[[which.i]] <- pick$qtl
        qtl.x <- gsub("Chr\\.", "X.", qtl.list[[which.i]])
        tmp <- strsplit(qtl.x[1], split="\\.")
        qtl[which.i] <- last.qtl <- paste(tmp[[1]][1],tmp[[1]][2], tmp[[1]][3], sep=".")
        oint[[which.i]] <- pick$oint
        blups[[which.i]] <- pick$blups
        if (is.null(add.qtl$xtra)) {
            phenoData[qtl.x] <- mD$asdata[qtl.list[[which.i]]] * 100
        }
        else {
            tmp.1 <- cbind.data.frame(geneticData[, merge.by],
                geneticData[, qtl.list[[which.i]]])
            names(tmp.1) <- c(merge.by, qtl.x)
            phenoData <- cbind.data.frame(ord = 1:nrow(phenoData),
                phenoData)
            phenoData <- merge(phenoData, tmp.1, by = merge.by,
                all.x = TRUE, all.y = FALSE)
            phenoData <- phenoData[order(phenoData$ord), ]
            phenoData <- phenoData[, -2]
        }
        gbind <- (2:dim(geneticData)[2])[!as.logical(state)]
        gD <- geneticData[, -gbind]
        mD <- mpmergeData(phenoData, gD, merge.by)
        cnt <- mD$cnt
        asdata <- mD$asdata
        qtls.cnt <- grep(paste(last.qtl, ".", sep=''), names(asdata), fixed=TRUE)
        which.i <- which.i + 1
        update <- TRUE
    }
    sigma2 <- add.qtl$sigma2
    if (names(add.qtl$gammas.con[length(add.qtl$gammas.con)]) ==
        "Fixed")
        sigma2 <- 1
    if (length(qtl)) {
        add.qtl$QTL$type <- gen.type
        add.qtl$QTL$qtl <- qtl
        add.qtl$QTL$qtl.list <- qtl.list
        add.qtl$QTL$effects <- cl[[which.i - 1]]
        add.qtl$QTL$veffects <- vl[[which.i - 1]]
        add.qtl$QTL$nfounders <- nfounders
        add.qtl$QTL$diag$dmat <- dmat
        add.qtl$QTL$diag$vl <- vl
        add.qtl$QTL$diag$cl <- cl
        add.qtl$QTL$diag$blups <- blups
        add.qtl$QTL$diag$oint <- oint
        add.qtl$QTL$state <- state
    }
    add.qtl <- mp.envDestruct(add.qtl, keep = c("asdata", "qtl.list",
        "ftrace", "data.name"))
    class(add.qtl) <- c("mpwgaim", "wgaim", "asreml")
    if(is.null(data.name)) {
        data.name <- paste(as.character(add.qtl$call$fixed[2]),
                           "data", sep = ".")
    }
    add.qtl$call$data <- as.name(data.name)
    cat(" Saving final data frame as ", data.name, " in working directory\n")
    assign(data.name, asdata, envir = parent.frame())
    add.qtl
}

