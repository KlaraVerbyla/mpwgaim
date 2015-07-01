##  Method remlrt.default
##
##' Default method for remlrt
##'
##' Default method which returns an error message
##' @title Default method for remlrt
##' @param object An object of class \code{mvmpwgaim}
##' @param \ldots Optional additional arguments.
##' @return Error message.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##'
remlrt.default <- function(object, ...) {
  stop("Currently the only supported method is \"mvwgaim\"")
}
