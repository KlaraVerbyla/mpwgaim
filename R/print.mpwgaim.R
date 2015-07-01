##  Method print.mpwgaim
##
##' Print method for mpwgaim objects
##'
##' Provides a simple printing of the selected QTL after using
##' 'mpwgaim'.
##' @title Print method for 'mpwgaim'
##' @param x an object of class 'mpwgaim'
##' @param intervalObj an object of class 'interval'
##' @param \ldots Additional optional objects
##' @author Julian Taylor (julian.taylor at adelaide.edu.au) and Ari
##' Verbyla (ari.verbyla at csiro.au)
##' @export
##'
print.mpwgaim <- function(x, intervalObj, ...)
{
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (is.null(x$QTL$qtl))
        cat("There are no significant putative QTL's\n")
    else {
        qtlm <- mpgetQTL(x, intervalObj)
        for (z in 1:nrow(qtlm)) {
            int <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
            if (x$QTL$type == "interval")
                cat("\nPutative QTL found on the interval", int,
                  "\nLeft-hand marker is", qtlm[z, 3], "\nRight-hand marker is",
                  qtlm[z, 5], "\n")
            else cat("\nPutative QTL found close to marker",
                int, "\nMarker is", qtlm[z, 3], "\n")
        }
    }
}
