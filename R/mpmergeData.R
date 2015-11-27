## Function mpmergeData
##
## Internal function for merging phenotypic and genetic data.
##
## Phenotypic data and the manipulated genetic data from an interval
## object are merged using the common indicator specified by 'by'.
## @title Merge phenotypic and genetic data
## @param phenoData Phenotypic data
## @param geneticData genetic data
## @param by the indicator variable on which to merge
## @return A list with the merged data frame ('asdata') and auxiliary
## information required for 'mpwgaim', including the columns of
## genetic data in 'asdata', columns for the selected QTL, and
## components relevant for dimension reduction.
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mpmergeData <- function (phenoData, geneticData, by) {
    int.cnt <- 2:dim(geneticData)[2]
    p <- length(int.cnt)
    whg <- geneticData[[by]] %in% unique(phenoData[[by]])
    ids <- as.character(geneticData[whg, by])
    q <- length(ids)
    phenoData <- cbind(ord = 1:nrow(phenoData), phenoData)
    geneticData[, int.cnt] <- geneticData[, int.cnt]/100
    if (p > q) {
##        require(MASS)
        mats <- geneticData[whg, int.cnt]
        tmat <- t(mats)
        xsvd <- svd(crossprod(tmat))
        xsvd.half <- t(xsvd$v %*% (t(xsvd$u) * sqrt(xsvd$d)))
        xsvd.inv <- ginv(xsvd.half)
        xtra <- tmat %*% xsvd.inv
        xsvd.df <- as.data.frame(xsvd.half)
        names(xsvd.df) <- paste("Tint.", 1:q, sep = "")
        xsvd.df[[by]] <- ids
        asdata <- merge(phenoData, xsvd.df, all.x = TRUE, by = by)
        asdata <- asdata[order(asdata$ord), ]
        asdata <- asdata[, -2]
        cnt <- grep("Tint\\.", names(asdata))
    }
    else {
        xtra <- NULL
        asdata <- merge(phenoData, geneticData, by.x = by, by.y = by,
            all.x = TRUE, all.y = FALSE)
        asdata <- asdata[order(asdata$ord), ]
        asdata <- asdata[, -2]
        dnams <- names(asdata)
        cnt <- grep("Chr\\.", dnams)
    }
    invisible(list(asdata = asdata, cnt = cnt,
        p = p, q = q, xtra = xtra))
}
