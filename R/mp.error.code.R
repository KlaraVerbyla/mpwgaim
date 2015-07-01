## Function mp.error.code
##
## Internal function to generate error message
##
## Generates an error message when the number of attempts to overcome
## convergence problems exceeds that set in \code{mvmpwgaim.asreml]
## or \code{remlrt}.
## @title printing of error messages on convergence failure
## @param ec Indicator of the convergence problem.
## @return NULL.  Printed message is the only output.
## @author Julian Taylor (julian.taylor at adelaide.edu.au)
##
mp.error.code <- function(ec = NULL){
  if(is.null(ec))
    stop("mpwgaim has been halted!")
  switch(ec,
         unstable = message("Error message:\nLikelihood not converging
or one or more parameters are unstable,
\ntry increasing the number of attempts or
change the model.  For diagnostic purposes the current fixed
\nQTL model is returned and full details of the current working random
QTL model can be found in ..$QTL"),
         converge = message("Warning message:\nParameter(s) not converging:
one or more parameters are unstable. Continuing with QTL analysis ....\n"))
}
