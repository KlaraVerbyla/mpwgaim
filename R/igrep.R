##  mpwgaim internal functions
##
## These functions are used internally in mpwgaim.
##
## igrep(patterns, x, ...)
## lgrep(patterns, x, ...)
## @title Internal mpwgaim functions.
## @param patterns list of character vectors.
## @param x object for which patterns are sought.
## @param \ldots Additional parameters.
## @return The index of matches of patterns in x.
## @author Julian Taylor <julian.taylor@adelaide.edu.au>
##
igrep <- function(patterns, x, ...){
    if(length(patterns) == 1)
        grep(patterns[[1]], x, ...)
    else {
        xt <- x
        for(i in 1:length(patterns)){
            ind <- grep(patterns[[i]], x, ...)
            x <- x[ind]
        }
        pmatch(x, xt)
    }
}
