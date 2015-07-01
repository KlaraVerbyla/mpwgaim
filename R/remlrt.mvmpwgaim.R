## Method remlrt.mvmpwgaim
##
##' Function to examine the significance of QTL by environment
##' interaction for each QTL
##'
##' Each QTL by environment term is removed in turn to examine if the
##' interaction is significant using a residual likelihood ratio test.
##' Those interactions that are not significant are removed from the
##' final model.  The original model used must have
##' \code{object$QTL$main} set to \code{TRUE}.
##' @title Residual likelihood ratio tests for QTL by environment interaction
##' @param object The result of using \code{mvmpwgaim.asreml} to find QTL.
##' @param data the data frame contained the data saved after using
##' \code{mvmpwgaim.asreml}.
##' @param attempts Maximum number of updates to carry out if there
##' are problems with convergence of the log-likelihood or parameter
##' estimates in fitting models using \code{asreml}.  Default is 5.
##' @param TypeI The type I error rate required.  Default is 0.05.
##' @param \ldots Other optional arguments.
##' @return A model with non-significant QTL by environment effects
##' removed and a table of tests in \code{object$QTL$REMLRT}.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
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
remlrt.mvmpwgaim <- function(object, data, attempts = 5, TypeI=0.05, ...) {
  ############### If main effects == TRUE and method = "random" , need to
  #####   (a) conduct REML tests for interactions
  #####   (b) remove interaction effects that are not significant
  #####   (c) remove main effects for those terms with interactions
  #####   (d) refit this reduced model
  ###############
#  if(object$QTL$method != "random") stop('REMLR tests only appropriate if QTL effects selected are random')
  if(!object$QTL$main) stop('Main effects not fitted in the analysis: refit with main.effects = TRUE')
  grp1 <- grep("X\\.", object$factor.names)
  which1 <- object$factor.names[grp1]
  grp2 <- grep(object$QTL$trait, which1)
  which2 <- which1[grp2]
  which2.str <- strsplit(which2, split=":")
  which3 <- c()
  #for (ii in 1: length(which2.str)) {
  #    which.ii <- grep("X\\.", which2.str[[ii]])
  #    which.jj <- 2
  #    if (which.ii == 2) which.jj <- 1
  #    #if(which.ii == 1) which3[ii] <- paste(paste("grp(", which2.str[[ii]][which.ii], ")", sep=""), ":", which2.str[[ii]][which.jj], sep="")
  #    #else which3[ii] <- paste(which2.str[[ii]][which.jj], ":", paste("grp(", which2.str[[ii]][which.ii], ")", sep=""), sep="")
  #    if(which.ii == 1) which3[ii] <- which2[ii] 
  #    else which3[ii] <- which2[[ii]]
  #}
  which3 = which2
  n.trait <- length(object$QTL$trait.levels)
  remlrt <- pvalue <- c()
  for (ii in 1:length(which2)) {
#    cat("\n Removing interaction: ", which2[ii], "in term ", which3[ii], "\n")
    cat("\n Removing interaction term ", which3[ii], "\n")
    cat("=======================\n")
    ran.form <- as.formula(paste("~ . - ", which3[ii], sep=""))
    loglik0 <- updateMPWgaim(object, data, attempts = attempts, random. = ran.form, ...)$loglik
    remlrt[ii] <- 2*(object$loglik - loglik0)
    pvalue[ii] <- (1-pchisq(remlrt[ii], df=1))/2
    if(remlrt[ii] < 0) {
      remlrt[ii] <- 0
      pvalue[ii] <- 1
    }
  }
  remlrt.df <- data.frame(Term = which2, Statistic = remlrt, Pvalue = pvalue)
  cat(" REMLRT statistics\n")
  cat("===================\n")
  qtl.type <- rep("interaction", length(which3))
  ind <- pvalue > TypeI
  which3 <- which3[ind]
  qtl.type[ind] <- "main"
  new.object <- object
  if(sum(ind) > 0) {
      tmp.form <- deparse(object$call$random[[2]], width.cutoff = 500)
      form.effects <- character(0)
      for (ii in 1:length(tmp.form)) {
          form.effects <- c(form.effects, strsplit(tmp.form[ii], split = " \\+ ")[[1]])
      }
    which3.out <- lgrep(as.list(which3), form.effects, fixed=TRUE)
    ran.effects <- form.effects[-which3.out]
    ran.term <- paste(ran.effects, collapse='+')
    ran.form <- as.formula(paste("~", ran.term, sep=""))
    cat("\n Updating model removing non-significant interactions .... \n")
    new.object <- updateMPWgaim(object, data, attempts = attempts, random. = ran.form, ...)
    new.object$QTL <- object$QTL
    new.object$QTL$qtl.type <- qtl.type
    list.coefs <- new.object$coefficients$random
    zind <- grep("X\\.", names(list.coefs))
    which.i <- length(object$QTL$diag$cl)+1
    new.object$QTL$effects <- new.object$QTL$diag$cl[[which.i]] <- list.coefs[zind]
    new.object$QTL$veffects <- new.object$QTL$diag$vl[[which.i]] <- new.object$vcoeff$random[zind]
    new.object$call$data <-  as.name(deparse(substitute(data)))
    class(new.object) <- class(object)
  }
  new.object$QTL$REMLRT <- remlrt.df
  new.object
}
