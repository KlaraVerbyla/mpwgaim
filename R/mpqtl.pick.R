## Function mpqtl.pick
##
##' Function for selecting a putative QTL in the forward selection
##' approach of mpwgaim.
##'
##' This function is usually only for internal use but could be used
##' to provide a manual rather than automated QTL analysis for
##' multi-parent populations.
##' @title QTL selection for 'mpwgaim'.
##' @param asr current QTL model.
##' @param intervalObj interval object constructed using
##' 'mpcross2int'.
##' @param asdata data frame containing the current phenotypic and
##' genetic information.
##' @param gen.type 'interval' or 'marker' depending on the type of
##' analysis.
##' @param state zero-one vector indicating the current QTL (ones).
##' @param verboseLev controls the printed output from this function;
##' any value greater than zero will flag that printing of outlier
##' statistics is to be carried out.
##' @return A list object with updated state vector, the qtl selected,
##' the outlier statistics for each interval or marker, and the blups
##' for each interval or marker.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##'
mpqtl.pick <- function(asr, intervalObj, asdata, gen.type, state, verboseLev)
{
    asr$call$data <- quote(asdata)
    sigma2 <- asr$sigma2
    if (names(asr$gammas.con[length(asr$gammas.con)]) == "Fixed")
        sigma2 <- 1
    if (!is.null(xtra <- asr$xtra)) {
        whq <- unlist(eval(asr$call$group$ints))
        n.int <- length(whq)
        cat(" Predict step for outlier statistics \n")
        cat("=====================================\n")
        qlev <- diag(n.int)
        qlist <- list(as.vector(qlev))
        names(qlist) <- "ints"
        pv <- predict(asr, classify = 'grp("ints")', only = 'grp("ints")',
            levels = qlist, vcov = TRUE, maxiter = 1, data = asdata)
        gamma <- pv$gammas[grep("ints", names(asr$gammas))]
        atilde <- pv$predictions$pvals[, "predicted.value"]
        atilde <- xtra %*% atilde
        pev <- pv$predictions$vcov
        vatilde <- sigma2 * gamma * diag(n.int) - pev
        vatilde <- apply(xtra, 1, function(el, vatilde) sum(el *
          (vatilde %*% el)), vatilde = vatilde)
        gnams <- names(state)[as.logical(state)]
    }
    else {
        gamma <- asr$gammas[grep("ints", names(asr$gammas))]
        rnams <- names(asr$coefficients$random)
        grp <- grep("ints", rnams)
        atilde <- asr$coefficients$random[grp]
        pevar <- sigma2 * asr$vcoeff$random[grp]
        vatilde <- sigma2 * gamma - pevar
        gnams <- substring(rnams[grp], first = 6)
    }
    names(atilde) <- gnams
    names(vatilde) <- gnams
    ntj2 <- ifelse(!is.na(atilde^2/vatilde), atilde^2/vatilde,
        0)
    names(ntj2) <- gnams
    split.gnams <- strsplit(gnams, split='\\.')
    mgnams <- unique(unlist(lapply(split.gnams,
                                   function(x) paste(x[1], x[2], x[3], sep='.'))))
    nisumtj2 <- c()
    noint <- length(mgnams)
    for (ni in 1:noint) {
        wh.m <- gnams[grep(mgnams[ni], gnams)]
        nisumtj2[ni] <- sum(atilde[wh.m]^2)/sum(vatilde[wh.m])
    }
    names(nisumtj2) <- mgnams
    nisumtj2[is.na(nisumtj2)] <- 0
    ##
    qtlint <- names(nisumtj2)[nisumtj2 == max(nisumtj2)]
    qsp <- unlist(strsplit(qtlint, split = "\\."))
    wint <- as.numeric(qsp[3])
    wchr <- qsp[2]
    oint <- c(nisumtj2)
    blups <- state
    qtl.interval <- paste("Chr", wchr, wint, sep='.')
    qtl.interval <- paste(qtl.interval, ".", sep='')
    qtl <- nint <- gnams[grep(qtl.interval, gnams, fixed = TRUE)]
    tint <- state[unique(unlist(lapply(split.gnams,
                                       function(x) paste(x[1], x[2], x[3], 1, sep='.'))))]
    tint[as.logical(tint)] <- oint
    oint <- tint
    blups[as.logical(state)] <- atilde/sqrt(vatilde)
    if (verboseLev > 0) {
        cgen <- "Interval"
        if (gen.type == "marker")
            cgen <- "Marker"
        cat(cgen, "outlier statistics \n")
        cat("=============================================== \n")
        for (i in 1:length(nisumtj2)) cat(cgen, names(nisumtj2)[i], "Outlier Statistic ",
            nisumtj2[i], "\n")
        cat("=============================================== \n\n")
    }
    state[nint] <- 0
    message("Found QTL on chromosome ", wchr, " ", gen.type,
        " ", wint)
    list(state = state, qtl = qtl, oint = oint, blups = blups)
}
