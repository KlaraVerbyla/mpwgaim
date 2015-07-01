## Function getmvQTL
##
## Internal function to obtain information for \code{summary.mvmpwgaim}
##
## This function extracts information from the map for selected QTL
## found using \code{mvmpwgaim}.
## @title Obtain map information for selected QTL
## @param object  An object of class \code{mvmpwgaim}.
## @param intervalObj  An interval object from using \code{mpcross2int}.
## @return A matrix of basic location information for each selected QTL.
## @author Julian Taylor (julian.taylor at adelaide.edu.au) and Ari
## Verbyla (ari.verbyla at csiro.au)
##
getmvQTL <- function(object, intervalObj)
{
  tmp.split <- strsplit(unlist(object$QTL$qtl), split="\\.")
  wchr <- unlist(lapply(tmp.split, function(el) el[2]))
  wint <- as.numeric(unlist(lapply(tmp.split, function(el) el[3])))
  qtlm <- matrix(ncol = 6, nrow = length(wchr))
  for (i in 1:length(wchr)) {
    lhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
    qtlm[i, 1:4] <- c(wchr[i], wint[i], names(lhmark),
                      round(lhmark,2))
    if (object$QTL$type == "interval") {
      if (length(intervalObj$geno[[wchr[i]]]$map) > 1)
        rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i] + 1]
      else rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
      qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
    }
    else qtlm <- qtlm[, -c(5:6)]
  }
  qtlm
}

