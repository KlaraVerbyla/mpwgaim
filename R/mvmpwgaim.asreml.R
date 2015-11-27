##  Method mvmpwgaim.asreml
##
##' Multivariate (multi-trait or multi-environment) QTL analysis for
##' multi-parent populations
##'
##' This function allows either a multi-trait or multi-environment QTL
##' analysis to be carried out.  Note that the current code does not
##' allow multiple traits and multiple environments in the same
##' analysis.  Two base models are required.  In the first, the
##' polygenic effect MUST have a diagonal structure for the trait or
##' environment variable, while in the second the polygenic effect
##' MUST have a factor analytic structure for the trait or environment
##' variable.  The diagonal form is used to establish if there is
##' sufficient evidence to select a putative QTL while the factor
##' analytic form is used both for selection and computation of impact
##' of the QTL through probability measures, log probability measures
##' and percentage variance accounted for by each QTL.
##'
##' For multi-environment analysis the logical argument main.effects
##' should be set to \code{TRUE}.  This ensures the possibility of
##' testing for environment by QTL interaction (see \code{remlrt}).
##'
##' Analysis using this function is both time consuming and
##' computationally demanding.  To facilitate analyses that take along
##' time, there is a logical variable that can be set, dorestart.  If
##' this argument is set to \code{TRUE}, a list of structures is saved
##' to file after each selection of a QTL.  If the analysis terminates
##' prematurely, this saved file can be used as input in a subsequent
##' run of the analysis through the restart argument.  The name of the
##' file saved is \code{"restart.RData"} and this should be the
##' argument in a subsequent call to \code{mvmpwgaim}.
##'
##' @title QTL analysis for multivariate multi-parent data
##' @param baseDiag A base model without marker based effects in
##' which all genetic effects associated with the multivariate nature
##' of the data are modelled using the 'diag' structure of 'asreml'
##' (see Details and Examples).
##' @param baseModel A base model without marker based effects in
##' which all genetic effects associated with the multivariate nature
##' of the data are modelled using the 'fa' structure of 'asreml'
##' (see Details and Examples).
##' @param phenoData Phenotypic data that includes all experimental
##' design and genetic line information that enables the base models
##' to be fitted.
##' @param intervalObj  An object of class 'interval'.
##' @param Trait A character vector that specifies the name of the
##' multivariate trait or environment variable in the phenotypic data.
##' @param merge.by A character vector that specifies the genetic
##' line variable that is present in both the phenotypic data and the
##' 'interval' object that is used in merging the two together.
##' @param gen.type A character vector indicating the type of
##' analysis to be undertaken.  Possible values are 'interval' (the
##' default) or 'marker'.
##' @param n.fa A number specifying the number of factors to be
##' included in the working model using founder probabilities.  The
##' default is 1 and the maximum is 2.  Note that the number can be
##' restricted by the number of traits or environments in the analysis.
##' @param TypeI The type I error rate for the forward selection of
##' putative QTL.  The default is 0.05.
##' @param attempts An integer representing the number of attempts at
##' convergence for the QTL model. The default is 5.
##' @param data.name A character vector specifying the name of the data
##' file to which the final data frame is to be saved.  If
##' \code{NULL}, the data frame is saved to the name of the response
##' dot data.
##' @param trace Logical for writing out debugging information.
##' Default is FALSE.
##' @param verboseLev An integer greater than 0 will result in
##' printing of outlier statistics in the QTL selection process.
##' @param main.effects Logical.  Should be \code{TRUE} for
##' multi-environment QTL analysis.  Otherwise the default is
##' \code{FALSE}.
##' @param restart A list with objects that allow the QTL analysis to
##' restart from the last successful QTL selection (see Details).
##' @param dorestart Logical.  Default is \code{FALSE}.  If
##' \code{TRUE}, the components required for a restart are saved after
##' each selection of a putative QTL.  This allows the process to be
##' restarted in the event that computational or time limits are
##' reached in the current analysis.
##' @param ... Additional objects that can be passed to functions, for
##' example 'asreml'.
##' @return A list object that has class 'mvmpwgaim' as well as
##' 'asreml'.  There is an additional list, \code{QTL}, which has
##' many components with information on the multi-parent population,
##' arguments set, QTL and diagnostics.
##' @author Ari Verbyla (ari.verbyla at csiro.au) and Klara Verbyla
##' (klara.verbyla at csiro.au)
##' @references Verbyla, A. P., Cavanagh, C. C. and Verbyla,
##' K. L. (2014). Whole genome analysis of multi-environment or
##' multi-trait QTL in MAGIC.  Genes, Genomes, Genetics.
##' @export
mvmpwgaim.asreml <-
  function (baseDiag, baseModel, phenoData, intervalObj, Trait = NULL, merge.by = NULL,
            gen.type = "interval", n.fa = 1, TypeI = 0.05, attempts = 10,
            data.name = NULL, trace = FALSE, verboseLev = 0, main.effects = FALSE,
            restart=NULL, dorestart=FALSE, ...)
{
    if(is.null(restart)) {
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
      if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
      if (is.null(other <- intervalObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
             "\".")
      if (is.null(phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
             "\".")
      if(n.fa > 2) stop("Number of factors for QTL model cannot be greater than 2.")
      mby <- pmatch(as.character(other), as.character(phenoData[,
                                                                merge.by]))
      if (all(is.na(mby)))
        stop("Names in Genotypic \"", merge.by, "\" column do not match any names in Phenotypic \"",
             merge.by, "\" column.")
      if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
      }
      if(intervalObj$gen.type == "mpPost") gen.type = "marker"
      if (gen.type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$intval)
      else gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
      if(intervalObj$gen.type != "mpPost") {
        if(any(lapply(gdat, function(x) dim(x)[2] %% nfounders) != 0))
          stop(" Genetic data is not consistent with the number of founders")
        nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
        lint <- unlist(lapply(nint, length))
        gnams <- paste("Chr", rep(names(intervalObj$geno), times = lint),
                       unlist(nint), sep = ".")
        eindex <- rep(1:nfounders,length(gnams))
        gnams <- paste(rep(gnams,each=nfounders), eindex, sep=".")
      }
      if(intervalObj$gen.type == "mpPost") {
        gnams <- unlist(lapply(intervalObj$geno, function(el) dimnames(el$imputed.data)[[2]]))
      }
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
        stop("Phenotypic data is not in same order as baseModel data.\n Try reordering phenoData apppropriately\n")
      n.trait <- length(levels(phenoData[, Trait]))
      if(n.trait > 2) {
        n.par.fa <- (n.fa+1)*n.trait - n.fa*(n.fa-1)/2
        n.par.us <- n.trait*(n.trait+1)/2
        if(n.par.fa > n.par.us) stop('n.fa set too high: reset and try again\n')
      }
      int.cnt <- 2:dim(geneticData)[2]
      mD <- mpmergeData(phenoData, geneticData, merge.by)
      cnt <- mD$cnt
      asdata <- mD$asdata
      state <- rep(1, length(int.cnt))
      names(state) <- names(geneticData[, int.cnt])
      if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
      }
      if (!baseDiag$converge) {
        cat("Warning: Base diagonal model has not converged. Updating base diagonal model\n")
        baseDiag <- update(baseDiag)
      }
      message("\nInvestigating the presence of a QTL......")
      cat("\n Diagonal Random Effects QTL plus Diagonal Polygenic Model\n")
      cat("===========================================================\n")
      qtlDiag <- baseDiag
      qtlDiag$call$group$ints <- cnt
      diag.term <- paste("diag(", Trait, "):grp('ints')",sep = "")
      diag.form <-  as.formula(paste("~ ",diag.term, " + .",sep = ""))
      qtlDiag <- updateMPWgaim(qtlDiag, asdata, attempts, random. = diag.form, ...)
      LRT <- 2*(qtlDiag$loglik - baseDiag$loglik)
      pvalue <- 1 - pchisq.mixture(LRT, ntrait = n.trait)
      if(pvalue < 0) pvalue <- 0
      cat("\nLikelihood Ratio Test Statistic: ", LRT, ", P-value: ", pvalue,"\n")
      dmat <- data.frame(L0 = baseDiag$loglik, L1 = qtlDiag$loglik, Statistic = LRT, Pvalue = pvalue)
      if(pvalue > TypeI) {
        cat("No QTL to be selected - no significant Trait QTL variation\n")
        return()
      }
      message("\nDeveloping the QTL model......")
      if (n.trait > 2) {
        cat("\n Preliminary fit for starting values\n")
        cat("=====================================\n")
        diagfa.qtl <- baseModel
        diagfa.qtl$call$group$ints <- cnt
        diagfa.qtl <- updateMPWgaim(diagfa.qtl, asdata, attempts, random. = diag.form, ...)

        cat("\n Factor analytic (1 factor) Random Effects QTL Model\n")
        cat("=====================================================\n")
        init.psi <- diagfa.qtl$gammas[grep('ints', names(diagfa.qtl$gammas))]*0.1
        init.lambda <- sqrt(diagfa.qtl$gammas[grep('ints', names(diagfa.qtl$gammas))])*0.9
        fa.init <- c(init.psi, init.lambda)
        names(fa.init) <- c(rep('P', n.trait), rep('U', n.trait))
        assign("fa.init", fa.init, envir = .GlobalEnv)
        add.qtl <- diagfa.qtl
        fa1 <- paste("fa(", Trait, ",1)", sep = "")
        fa1.term <- paste(fa1, ":grp('ints')", sep = "")
        new.random <- as.formula(paste("~ ", fa1.term, "+ . -", diag.term, sep = ""))
        add.qtl$call$random <- update.formula(add.qtl$call$random, new.random)
        add.qtl$random <- update.formula(add.qtl$random, new.random)
        names(add.qtl$G.param)[1] <- fa1
        add.qtl$G.param[[1]][[1]]$initial <- fa.init
        fa1.nam1 <-  paste(paste(fa1, ":", merge.by, "!", Trait, sep=""), add.qtl$G.param[[1]][[1]]$levels, sep=".")
        fa1.nam2 <- paste(fa1.nam1, rep(c("var", "fa1"), each = length(add.qtl$G.param[[1]][[1]]$levels)), sep=".")
        names(add.qtl$G.param[[1]][[1]]$initial) <- fa1.nam2
        names(add.qtl$G.param[[1]][[2]]$initial) <- gsub(Trait, paste("fa(", Trait," 1)", sep=""), names(add.qtl$G.param[[1]][[2]]$initial))
        add.qtl$G.param[[1]][[1]]$con <- rep(c("P", "U"), each = length(add.qtl$G.param[[1]][[1]]$levels))
        add.qtl$G.param[[1]][[1]]$model <- "fa"
        add.qtl$G.param[[1]][[1]]$levels <- c(add.qtl$G.param[[1]][[1]]$levels, "Comp1")
        add.qtl <- updateMPWgaim(add.qtl, asdata, attempts, ...)
        if (n.fa == 2) {
          cat("\n Factor analytic (2 factors) Random Effects QTL Model\n")
          cat("======================================================\n")
          init.psi <- add.qtl$gammas[grep('ints', names(add.qtl$gammas))][1:n.trait]*0.1
          fa1.gammas <- add.qtl$gammas[igrep(list('fa1','ints'), names(add.qtl$gammas))]
          fa2.init <- c(init.psi, fa1.gammas*0.9, 0, fa1.gammas[2:n.trait]*0.1)
          fa2 <- paste("fa(", Trait, ",2)", sep = "")
          fa2.term <- paste(fa2, ":grp('ints')", sep = "")
          new.random <- as.formula(paste("~ ", fa2.term, "+ . -", fa1.term, sep = ""))
          add.qtl$random <- update.formula(add.qtl$random, new.random)
          names(add.qtl$G.param)[1] <- fa2
          add.qtl$G.param[[1]][[1]]$initial <- fa2.init
          fa2.nam1 <-  paste(paste(fa2, ":", merge.by, "!", Trait, sep=""), add.qtl$G.param[[1]][[1]]$levels[-length(add.qtl$G.param[[1]][[1]]$levels)], sep=".")
          fa2.nam2 <- paste(fa2.nam1, rep(c("var", "fa1", "fa2"), each = length(add.qtl$G.param[[1]][[1]]$levels))-1, sep=".")
          names(add.qtl$G.param[[1]][[1]]$initial) <- fa2.nam2
          names(add.qtl$G.param[[1]][[2]]$initial) <- gsub(paste("fa(", Trait," 1)", sep=""), paste("fa(", Trait," 2)", sep=""), names(add.qtl$G.param[[1]][[2]]$initial))
          add.qtl$G.param[[1]][[1]]$con <- c(add.qtl$G.param[[1]][[1]]$con, "F", rep("U", each = length(add.qtl$G.param[[1]][[1]]$levels)-1))
          add.qtl$G.param[[1]][[1]]$levels <- c(add.qtl$G.param[[1]][[1]]$levels, "Comp2")
          add.qtl <- updateMPWgaim(add.qtl, asdata, attempts, ...)
        }
      }
      if (n.trait == 2) {
          add.qtl <- baseModel
          add.qtl$call$group$ints <- cnt
          cat("\n Bivariate Correlation Random Effects QTL Model\n")
          cat("================================================\n")
          cor.term <-  paste("corh(", Trait, "):grp('ints')", sep = "")
          cor.form <- as.formula(paste("~ ",cor.term, " + .",sep = ""))
          add.qtl <- updateMPWgaim(add.qtl, asdata, attempts, random. = cor.form, ...)
          corg <- add.qtl$gammas[igrep(list("cor","ints"), names(add.qtl$gammas))]
      }
###################
      update <- FALSE
      qtl.int <- c()
      vl <- cl <- qtl <- ochr <- oint <- list()
      which.i <- 1
      ochr <- oint <- blups <- restart <- list()
    }
    else {
      update <- restart$update
      which.i <- restart$which.i
      qtl.int <- restart$qtl.int
      vl <- restart$vl
      cl <- restart$cl
      qtl <- restart$qtl
      oint <- restart$oint
      blups <- restart$blups
      nfounders <- restart$nfounders
      add.qtl <- restart$add.qtl
      baseDiag <- restart$baseDiag
      qtlDiag <- restart$qtlDiag
      LRT <- restart$LRT
      pvalue <- restart$pvalue
      dmat <- restart$dmat
      state <- restart$state
      phenoData <- restart$phenoData
      geneticData <- restart$geneticData
      Trait <- restart$Trait
      main.effects <- restart$main.effects
      merge.by <- restart$merge.by
      last.qtl <- restart$last.qtl
      cnt <- restart$cnt
      qtls.cnt <- restart$qtls.cnt
      asdata <- restart$asdata
      attempts <- restart$attempts
      n.trait <- restart$n.trait
      TypeI <- restart$TypeI
      intervalObj <- restart$intervalObj
      gbind <- (2:dim(geneticData)[2])[!as.logical(state)]
      gD <- geneticData[, -gbind]
      mD <- mpmergeData(phenoData, gD, merge.by)
    }
    message("Searching for QTLs......")
    repeat {
      if (update) {
        add.qtl$call$group[[last.qtl]] <- baseDiag$call$group[[last.qtl]] <- qtlDiag$call$group[[last.qtl]] <- qtls.cnt
        qtlDiag$call$group$ints <- add.qtl$call$group$ints <- cnt
        if(main.effects) {
          new.term <- paste(Trait, ":grp(",last.qtl,")", sep='')
          new.term <- paste(paste("grp(",last.qtl,")", sep=""),  new.term, sep = "+")
        }
        else {
          new.term <- paste("diag(", Trait, "):grp(",last.qtl,")", sep='')
        }
        new.form <- as.formula(paste("~ . +", new.term, sep = ""))
        cat("\nDiagonal QTL Model Iteration (", which.i, "):\n")
        cat("==================================\n")
        baseDiag <- updateMPWgaim(baseDiag, asdata, attempts, random. = new.form, ...)
        cat("\nDiagonal QTL and Interval/Marker Model Iteration (", which.i, "):\n")
        cat("======================================================\n")
        qtlDiag <- updateMPWgaim(qtlDiag, asdata, attempts, random. = new.form, ...)
        LRT[which.i] <- 2*(qtlDiag$loglik - baseDiag$loglik)
        pvalue[which.i] <- 1-pchisq.mixture(LRT[which.i], ntrait=n.trait)
        if(pvalue[which.i] < 0) pvalue[which.i] <- 0
        dmat[which.i, ] <- c(baseDiag$loglik, qtlDiag$loglik, LRT[which.i],
                             pvalue[which.i])
        cat("\nLikelihood Ratio Test Statistic: ", LRT[which.i], ", P-value: ", pvalue[which.i],"\n")
        if(pvalue[which.i] < TypeI) {
          cat("\nQTL Selection Iteration (", which.i,"):\n")
          cat("============================\n")
        }
        else {
          cat("\n Fitting final model:\n")
          cat("======================\n")
        }

        #
        # Check that the initial values for fa model terms have non-zero specific variances
        #
        if(n.trait > 2) add.qtl <- mp.fa.modify(add.qtl, Trait, merge.by)
        else add.qtl <- mpch.modify(add.qtl, Trait, merge.by)
        add.qtl <- updateMPWgaim(add.qtl, asdata, attempts, random. = new.form, ...)
        #cat(" Random model\n")
        #print(add.qtl$call$random)
        #cat(" Variance parameter estimates\n")
        #print(add.qtl$gammas)
        list.coefs <- add.qtl$coefficients$random
        zind <- grep("X\\.", names(list.coefs))
        list.coefs <- list.coefs[zind]
        cl[[which.i - 1]] <- list.coefs
        vl[[which.i - 1]] <- add.qtl$vcoeff$random[zind]
      }
      if(pvalue[which.i] > TypeI) break
      if(mD$p > mD$q)
        add.qtl$xtra <- mD$xtra
      labels <- c(Trait, "ints")
      pick <- mvmp.pick(add.qtl, intervalObj, asdata, gen.type, state, labels, n.trait, n.fa, verboseLev)

      state <- pick$state
      qtl[[which.i]] <- pick$qtl
      qtl.x <- gsub("Chr\\.", "X.", qtl[[which.i]])
      tmp <- strsplit(qtl.x[1], split="\\.")
      qtl.int[which.i] <- last.qtl <- paste(tmp[[1]][1],tmp[[1]][2], tmp[[1]][3], sep=".")
      oint[[which.i]] <- pick$oint
      blups[[which.i]] <- pick$blups
      if (is.null(add.qtl$xtra)) {
        phenoData[qtl.x] <- mD$asdata[qtl[[which.i]]] * 100
      }
      else {
        tmp.1 <- cbind.data.frame(geneticData[, merge.by],
                                  geneticData[, qtl[[which.i]]])
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
      if(dorestart){
        restart$update <- update
        restart$which.i <- which.i
        restart$qtl.int <- qtl.int
        restart$vl <- vl
        restart$cl <- cl
        restart$qtl <- qtl
        restart$oint <- oint
        restart$blups <- blups
        restart$nfounders <- nfounders
        restart$add.qtl <- add.qtl
        restart$baseDiag <- baseDiag
        restart$qtlDiag <- qtlDiag
        restart$LRT <- LRT
        restart$pvalue <- pvalue
        restart$dmat <- dmat
        restart$state <- state
        restart$phenoData <- phenoData
        restart$geneticData <- geneticData
        restart$Trait <- Trait
        restart$main.effects <- main.effects
        restart$merge.by <- merge.by
        restart$last.qtl <- last.qtl
        restart$cnt <- cnt
        restart$qtls.cnt <- qtls.cnt
        restart$asdata <- asdata
        restart$attempts <- attempts
        restart$n.trait <- n.trait
        restart$TypeI <- TypeI
        restart$intervalObj <- intervalObj
        filename = paste0(as.character(add.qtl$call$fixed[2]),'restart.RData')
        cat(" Saving current structures to", filename, "\n")
        save(restart, file= filename)
      }
    }
    sigma2 <- add.qtl$sigma2
    if (names(add.qtl$gammas.con[length(add.qtl$gammas.con)]) ==
          "Fixed")
      sigma2 <- 1
    qtl.list <- list()
    if (length(qtl)) {
      qtl.list$type <- gen.type
      qtl.list$qtl <- qtl.int
      qtl.list$qtl.list <- qtl
      qtl.list$nfounders <- nfounders
      qtl.list$diag$dmat <- dmat
      qtl.list$diag$vl <- vl
      qtl.list$diag$cl <- cl
      qtl.list$diag$blups <- blups
      qtl.list$diag$oint <- oint
      qtl.list$effects <- cl[[which.i - 1]]
      qtl.list$veffects <- vl[[which.i - 1]]
      if(main.effects) qtl.list$main <- TRUE
      else qtl.list$main <- FALSE
    }
    qtl.list$trait <- Trait
    qtl.list$trait.levels <- levels(asdata[,Trait])
    if(is.null(data.name)) {
      data.name <- paste(as.character(add.qtl$call$fixed[2]),
                         "data", sep = ".")
    }
    add.qtl$call$data <- as.name(data.name)
    cat(" Saving final data frame as ", data.name, " in working directory\n")
    assign(data.name, asdata, envir = parent.frame())
    final.model <- mp.envDestruct(add.qtl, keep = c("data.name", "qtl.list",
                                                 "ftrace"))
    final.model$QTL <- qtl.list
    class(final.model) <- c("mvmpwgaim", "mvwgaim", "asreml")
    final.model
}


## @examples
## \dontrun {
## ## Simulation of genetic and phenotypic data for analysis
## require(asreml)
## require(mpMap)
## require(qtl)
## ## Generate a linkage map
## map <- sim.map(len=rep(200,7), n.mar=rep(51,7), eq.spacing=TRUE, include.x=FALSE)
## ## Set up a pedigree  for a 4-way cross with 500 RILs and 6 generations of selfing
## sim.ped <- sim.mpped(4, 1, 500, 6, 1)
## ## Dummy QTL structure
## qtl.mat <- matrix(data=c(1, 142, 0.354, -0.354, -0.354,  0.354,
##                          2, 162, 0.354, -0.354, -0.354, -0.354,
##                          5, 78, 0.354, -0.354, -0.354, 0.354),
##                          nrow=3, ncol=6, byrow=TRUE)
## ## Simulate mpcross object
## sim.dat <- sim.mpcross(map=map, pedigree=sim.ped,
##                        qtl=qtl.mat, seed=5)
## mpSim <- maporder(sim.dat)
## ## Create interval objects (see mpcross2int)
## mpInterval <- mpcross2int(mpSim, gen.type="mpInterval")
## mpMarker <- mpcross2int(mpSim, gen.type="mpMarker")
## ##  Model matrix for computing the contribution of the QTL
## nqtl.col <- dim(sim.dat$qtlgeno$finals)[2]
## mmat <- matrix(nrow=500, ncol=4*nqtl.col/2)
## for (ii in 1:(nqtl.col/2)) {
## qtl.fac <- factor(sim.dat$qtlgeno$finals[,ii])
## mmat[,(4*ii-3):(4*ii)] <- model.matrix(~qtl.fac - 1)
## }
## ##  Effects for each environment
## ## Environment 1
## qtl.mat1 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##                           2, 162, 0.354, -0.354, -0.354, 0.354,
##                           5, 78, 0.354, -0.354, -0.354, 0.354),
##                           nrow=3, ncol=6, byrow=TRUE)
## qtl.sizes1 <- as.vector(t(qtl.mat1[,3:6]))
## qtl.effect1 <- mmat %*% qtl.sizes1
## ## Environment 2
## qtl.mat2 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##                           2, 162, 0.354, -0.354, -0.354, 0.354,
##                           5, 78, 0, 0, 0, 0),
##                           nrow=3, ncol=6, byrow=TRUE)
## qtl.sizes2 <- as.vector(t(qtl.mat2[,3:6]))
## qtl.effect2 <- mmat %*% qtl.sizes2
## ## Environment 3
## qtl.mat3 <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
##                           2, 162, -0.354, 0.354, 0.354, -0.354,
##                           5, 78, 0.354, -0.354, -0.354, 0.354),
##                           nrow=3, ncol=6, byrow=TRUE)
## qtl.sizes3 <- as.vector(t(qtl.mat3[,3:6]))
## qtl.effect3 <- mmat %*% qtl.sizes3
## ## Polygenic variance
## pvar_0.5
## ## Function to calculate approximate percentage variance for each QTL
## perc.var <- function(qtl.mat, poly.var) {
##    nfounders <- dim(qtl.mat)[2]-2
##    prob <- 1/nfounders
##    varq <- diag(rep(prob,nfounders)) - rep(prob,nfounders) %*% t(rep(prob,nfounders))
##    gvar <- apply(qtl.mat[, -c(1,2)], 1, function(el, varq) sum((el %*% varq) * el), varq)
##    totvar <- sum(gvar)+pvar
##    perc.var <- 100*gvar/totvar
##    round(perc.var,1)
## }
## ## Percentage variance for each QTL in each environment
## ## Environment 1
## percvar1 <- perc.var(qtl.mat1, pvar)
## percvar1
## ## Environment 2
## percvar2 <- perc.var(qtl.mat2, pvar)
## percvar2
## ## Environment 3
## percvar3 <- perc.var(qtl.mat3, pvar)
## percvar3
## ## Setup simulated data for analysis
## ntrait <- 3
## ngeno <- 500
## nrep <- 2
## Trait <- factor(rep(rep(1:ntrait, ngeno), nrep))
## gvar <- matrix(c(1.0,0.9,0.7,
##                  0.9,1.0,0.7,
##                  0.7,0.7,1.0), ncol=3)
## gchol <- t(chol(gvar))
## id <- factor(rep(rep(paste0("L", 1:ngeno),each=ntrait),nrep))
## pheno.data <- data.frame(y = 0, Trait=Trait, Rep=nrep, id=id)
## qtleffects <- as.matrix(0, ncol=1, nrow=ntrait*ngeno)
## qtleffects[rep(1:3,ngeno)==1] <- qtl.effect1
## qtleffects[rep(1:3,ngeno)==2] <- qtl.effect2
## qtleffects[rep(1:3,ngeno)==3] <- qtl.effect3
## set.seed(10)
## ee <- rnorm(ntrait*ngeno*nrep,0,1)
## geno <- rnorm(ntrait*ngeno,0,1)
## uu <- as.vector(gchol%*%matrix(geno, nrow=ntrait))
## pheno.data$y <- rep(c(9,10,12), ngeno*nrep)+
##     qtleffects + sqrt(1/2)*c(uu,uu) + ee
## ord <- order(pheno.data$Trait)
## pheno.data<-pheno.data[ord,]
## ## Multi-environment analysis
## ## Multi-environment preliminary models
## test.asr0 <- asreml(y ~ Trait - 1, random = ~ diag(Trait):id, data=pheno.data)
## test.asr1 <- asreml(y ~ Trait - 1, random = ~ fa(Trait,1):id, data=pheno.data)
## ##  MVMPWGAIM analysis
## test.qtl <- mvmpwgaim(test.asr0, test.asr1, pheno.data, mpInterval, merge.by = "id",
##                       Trait="Trait", verboseLev=0, main.effects = TRUE,
##                       data.name = "mvy.data", gen.type="interval", na.method.X='include',
##                       workspace=4e7, pworkspace=4e7)
## ## Check for Trait by QTL interactions and fit a reduced model if needed
## test.aqtl <- remlrt.mvmpwgaim(test.qtl, mvy.data, TypeI=0.05)
## ##  Summary for MVMPWGAIM analysis
## test.summ <- summary(test.aqtl, mpInterval)
## test.summ
## }
##
