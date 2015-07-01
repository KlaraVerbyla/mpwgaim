## Function mvmp.pick
##
##' Internal function for QTL selection in MAGIC
##'
##' This function is called from \code{mvmpwgaim.asreml} when
##' selection of a QTL is justified.
##' @title QTL selection for \code{mvmpwgaim.asreml}
##' @param asr fitted model containing marker-based information
##' @param intervalObj an interval object that contains the genetic
##' information necessary for selection of a putative QTL.
##' @param asdata data frame containing the phenotypic and genetic
##' data structure.
##' @param gen.type Either '"interval"' or '"marker"'
##' @param state Indicator variable of those markers or intervals not
##' yet selected.
##' @param labels Character vector with labels for the Trait or
##' environment factor and '"ints"', the internal label for the
##' genetic information contained in '"asdata"'.
##' @param n.trait Number of traits or environments
##' @param n.fa How many factors in the factor analytic model that
##' has been fitted in \code{asr}.
##' @param verboseLev Controls the printing of output from this
##' function.  Verbose printing if the value of this argument is
##' greater than 0.
##' @param debug Logical set at \code{FALSE}.  Used by the developers
##' to debug the code.
##' @return A list object is returned with the \code{state} vector
##' updated, the QTL selected, outlier statistics used in the
##' selection process and best linear unbiased predictors for the
##' potential effects in each interval or at each marker depending on
##' \code{gen.type}.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##'
mvmp.pick <- function(asr, intervalObj, asdata, gen.type, state, labels,
                      n.trait, n.fa, verboseLev, debug=FALSE) {
##    require(MASS) ## delete - ABZ
    asr$call$data <- quote(asdata)
    nfounders <- intervalObj$nfounders
    gnams <- names(state)[as.logical(state)]
    whq <- unlist(eval(asr$call$group[[labels[2]]]))  ##### problem
    if(debug) print(whq)
    n.int <- length(whq)
    if(debug) print(n.int)

    cat(" Predict step for outlier statistics \n")
    cat("=====================================\n")
#
#  use predict for random effects a and their pev
#
    qlev <- diag(n.int)
    qlist <- list(as.vector(qlev))
    names(qlist) <- labels[2]

    if(n.trait != 2) {
        pv <- predict(asr, classify = paste(labels[1], labels[2], sep=":"),
                      only = paste("fa(",labels[1], ", ", n.fa, "):", labels[2], sep=""),
                      levels=qlist, vcov=TRUE, maxiter = 1, data = asdata)
        if(debug) {
            print(pv$gammas)
            print(head(pv$predictions$pvals))
        }
        gams <- pv$gammas[igrep(list(labels[1], labels[2]), names(pv$gammas))]
        psi <- gams[1:n.trait]
        if(debug) print(psi)
        Lam <- matrix(gams[(n.trait + 1):((n.fa + 1) * n.trait)], ncol = n.fa)
        if(debug) print(Lam)
        sigma2 <- pv$sigma2
        if (names(pv$gammas.con[length(pv$gammas.con)]) == "Fixed")
            sigma2 <- 1
        Ga <- sigma2 * (Lam %*% t(Lam) + diag(psi))
        if(debug) {
            cat(" Ga\n")
            print(Ga)
        }
    }
    else {
        pv <- predict(asr, classify = paste(labels[1], labels[2], sep=':'),
                      only = paste0(labels[1], sep=':', 'grp(\"', labels[2], '\")'),
                      levels=qlist, vcov=TRUE, maxiter = 1, data=asdata)
        if(debug) print(head(pv$predictions$pvals))
        Ga <- summary(pv, nice = TRUE)$nice[[paste0(labels[1], sep=':', 'grp(\"', labels[2], '\")')]][['grp(\"ints\")']]
    ####### Note that corh does not give the correct Ga: need to modify nice structure
        Ga[1,2] <- Ga[2,1] <- Ga[1,2]*sqrt(Ga[1,1]*Ga[2,2])
        sigma2 <- pv$sigma2
        if (names(pv$gammas.con[length(pv$gammas.con)]) == "Fixed")
            sigma2 <- 1
        Ga <- sigma2*Ga
        if(debug) {
            print(sigma2)
            print(pv$gammas)
            cat(" Ga\n")
            print(Ga)
        }
    }
    ord <- order(pv$predictions$pvals[, labels[1]], pv$predictions$pvals[, labels[2]])
    pvals <- pv$predictions$pvals[ord,]
    if(debug) print(head(pvals))
    vcov <- pv$predictions$vcov[ord,ord]
    if(debug) print(sqrt(diag(vcov)))
    if(debug) print(xyplot(sqrt(diag(vcov)) ~ pvals[,4]))
    if (any(diag(vcov) < 0)) cat(" Negative variances\n")
    atilde <- pvals[, 'predicted.value']
    atilde[is.na(atilde)] <- 0
    matilde <- matrix(atilde, ncol=n.trait, byrow=FALSE)
    if(debug) {
        cat(" matilde\n")
        print(matilde)
    }
  #  covariance matrix
    Ginv <- ginv(Ga)
    if(debug) {
        cat(" Ginv\n")
        print(Ginv)
    }
    pev <- vcov
    pev[is.na(pev)] <- 0
    svd.pev <- svd(pev)
    if(debug) {
        cat(" pev\n")
        print(svd.pev$d)
    }
    vatilde <- kronecker(Ga, diag(n.int)) - pev
    if(debug) print(diag(vatilde))
    svd.vat <- svd(vatilde)
    if(debug) {
        cat(" svd vatilde\n")
        print(svd.vat$d)
    }
  #   outlier statistics
    if(!is.null(xtra <- asr$xtra)) {
        if(debug) print(dim(xtra))
        matilde <- xtra %*% matilde
    } else {
        xtra <- diag(dim(matilde)[1])
    }
    dimnames(matilde)[[1]] <- gnams
  #      numerators
    nums <- apply(matilde, 1, function(el, Ginv) sum(el * (Ginv %*% el)), Ginv = Ginv)
    if(debug) {
        cat(" nums\n")
        print(nums)
    }
                                        #      denominators
    denom <- apply(xtra, 1, function(el, Ginv, vatilde)
               {tmp1 <- kronecker(diag(n.trait), el)
                tmp2 <- t(tmp1) %*% vatilde %*% tmp1
                tmp3 <- sum(diag(Ginv %*% tmp2))
                if(debug) {
                    cat(" denom step\n")
                    print(tmp1)
                    print(tmp2)
                    if(tmp3 < 0) cat("This one\n")
                }
                tmp3
            },
                   Ginv = Ginv, vatilde = vatilde)
    if(debug) {
        cat(" denom\n")
        print(denom)
    }
    mat.nums <- matrix(nums, ncol=nfounders, byrow=TRUE)
    mat.denom <- matrix(denom, ncol=nfounders, byrow=TRUE)
    if(debug) print(mat.denom)
    tnums <- rowSums(mat.nums)
    tdenom <- rowSums(mat.denom)
    t2kj <- ifelse(!is.na(tnums/tdenom), tnums/tdenom, 0)
    split.gnams <- strsplit(gnams, split='\\.')
    tgnams <- unique(unlist(lapply(split.gnams,
                                   function(x) paste(x[1], x[2], x[3], sep='.'))))
    names(t2kj) <- tgnams
    tnint <- names(t2kj)[t2kj == max(t2kj)]
    qsp <- unlist(strsplit(tnint, split = "\\."))
    wint <- as.numeric(qsp[3])
    wchr <- qsp[2]
    if(verboseLev > 0){
        cgen <- "Interval"
        if (gen.type == "marker")
            cgen <- "Marker"
        cat(cgen, "outlier statistics \n")
        cat("=============================================== \n")
        for (i in 1:length(t2kj)) cat(cgen, names(t2kj)[i], "Outlier Statistic ",
                                      t2kj[i], "\n")
        cat("=============================================== \n\n")
    }
    oint <- c(t2kj)
    tint <- state[unique(unlist(lapply(split.gnams,
                                       function(x) paste(x[1], x[2], x[3], 1, sep='.'))))]
    tint[as.logical(tint)] <- oint
    oint <- tint
    split.onams <- strsplit(names(oint), split='\\.')
    tint.names <- unlist(lapply(split.onams,
                                function(x) paste(x[1], x[2], x[3], sep='.')))
    names(oint) <- tint.names
    qtl.interval <- paste(tnint, ".", sep='')
    qtl <- nint <- gnams[grep(qtl.interval, gnams, fixed = TRUE)]
    state[nint] <- 0
    if(!is.null(asr$xtra)) {
        blups.se <- c()
        ng <- dim(xtra)[2]
        for (ii in 1:n.trait) {
            vatilde.ii <- vatilde[((ii-1)*ng+1):(ii*ng),((ii-1)*ng+1):(ii*ng)]
            bse <- apply(xtra, 1,function(el, vatilde.ii) sum(el * (vatilde.ii %*% el)),
                                               vatilde.ii = vatilde.ii)
            if(any(bse < 0)) cat(" Negative standard errors!\n")
            blups.se <- c(blups.se, sqrt(bse))
        }
    }
    else blups.se <- sqrt(diag(vatilde))
    blups.se <- matrix(blups.se, ncol=n.trait, byrow=FALSE)
    if(debug) {
        cat(" blups.se\n")
        print(blups.se)
    }
    blups <- matilde/blups.se
    dimnames(blups)[[2]] <- levels(asdata[,labels[1]])
    message("Found QTL on chromosome ", wchr, " ", gen.type,
            " ", wint)
    list(state = state, qtl = qtl, oint = oint, blups = blups)
}
