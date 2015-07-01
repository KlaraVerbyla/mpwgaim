## Function print.summary.mpwgaim
##
##' Print method for objects of class 'summary.mpwgaim'
##'
##' Prints the summary information found using 'summary.mpwgaim'
##' @title Print table of QTL information
##' @param x an object of class 'summary.mpwgaim'
##' @param \ldots other arguments passed to print.summary.mpwgaim
##' (ignored)
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##'
print.summary.mpwgaim <- function(x,...) {
    print(x$summary)
    invisible()
}
