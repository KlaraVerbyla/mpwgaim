## Method function for mvmpwgaim
##
##' Use the appropriate method for using mvmpwgaim.  Currently only
##' the method 'asreml' is supported.
##' @title mvmpwgaim Method function
##' @param baseDiag Base model with a 'diag' structure for the
##' polygenic effects.
##' @param baseModel Base model with a factor analytic structure for
##' the polygenic effects.
##' @param \ldots Additional arguments.
##' @return AN object of class 'mvmpwgaim'
##' @author Ari verbyla (ari.verbyla at csiro.au)
##' @export
##'
mvmpwgaim.default <- function(baseDiag, baseModel, ...) {
  stop("Currently the only supported method is \"asreml\"")
}
