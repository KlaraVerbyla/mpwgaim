## Method mpwgaim.default
##
##' Internal mpwgaim function
##'
##' Default method for fitting mpwgaim.  The only available method is
##' 'asreml'.  This function simply issues a message and stops.
##' @title mpwgaim for method 'default'
##' @param baseModel a model object
##' @param \ldots Additional parameter
##'
mpwgaim.default <- function(baseModel, ...) {
    stop("Currently the only supported method is \"asreml\"")
}
