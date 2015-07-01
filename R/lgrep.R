##  mpwgaim internal functions
##
## These functions are used internally in mpwgaim.
##
## igrep(patterns, x, ...)
## lgrep(patterns, x, ...)
## @title Internal mpwgaim functions.
## @param patterns: list of character vectors.
## @param x: object for which patterns are sought.
## @param ... Additional parameters.
## @return The index of matches of patterns in x.
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
lgrep <- function(patterns, x, ...){
    if(length(patterns) == 1)
        grep(patterns[[1]], x, ...)
    else {
        ind <- list()
        for(i in 1:length(patterns)){
            ind[[i]] <- grep(patterns[[i]], x, ...)
        }
        sort(unlist(ind))
    }
}

