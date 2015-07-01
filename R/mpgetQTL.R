## Function mpgetQTL
##
##  Internal function to obtain basic QTL information
##
## Extracts information on the QTL used in summary.mpwgaim.
## @title  Get QTL information for 'summary.mpwgaim'
## @param object An object of class 'mpwgaim'
## @param intervalObj  An object of class 'interval'
## @return A matrix of marker information for each QTL.
## @author Julian Taylor (julian.taylor at adelaide.edu.au) and Ari
## Verbyla (ari.verbyla at csiro.au)
##
mpgetQTL <-  function(object, intervalObj)
  {
    spe <- strsplit(object$QTL$qtl, split = "\\.")
    wchr <- unlist(lapply(spe, function(el) el[2]))
    wint <- as.numeric(unlist(lapply(spe, function(el) el[3])))
    qtlm <- matrix(ncol = 6, nrow = length(wchr))
    for (i in 1:length(wchr)) {
        lhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
        qtlm[i, 1:4] <- c(wchr[i], wint[i], names(lhmark), round(lhmark, 2))
        if (object$QTL$type == "interval") {
            if (length(intervalObj$geno[[wchr[i]]]$map) > 1)
                rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i] +
                                                          1]
            else rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
            qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
        }
    }
    if (object$QTL$type == "marker") {
        qtlm <- qtlm[, -c(5:6), drop = FALSE]
    }
    qtlm
}
