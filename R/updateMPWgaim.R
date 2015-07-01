## Function updateMPWgaim
##
## Internal function for fitting models in 'mpwgaim.asreml'
##
## This function is used to fit models in 'mpwgaim' by allowing
## additional terms to be included.  It also traps errors and ensures
## convergence of the log-likelihood and parameter estimates is
## monitored.
## @title Model fitting in 'mpgaim.asreml'
## @param object 'asreml' object
## @param asdata data frame to be used in fitting the model
## @param attempts number of attempts to update to ensure both the
## log-likelihood and parameter estimates have converged.
## @param \ldots Optional additional objects to be passed to the
## update function.
## @return: Fitted model of class 'asreml'
## @author Julian Taylor (julian.taylor at adelaide.edu.au) and Ari
## Verbyla (ari.verbyla at csiro.au)
##
updateMPWgaim <- function(object, asdata, attempts, ...)
{
    object$call$data <- quote(asdata)
    class(object) <- "asreml"
    object <- update(object, ...)
    fwarn <- function(object) {
        mon <- object$monitor
        prgam <- mon[4:nrow(mon), (ncol(mon) - 2)]
        prgam[abs(prgam) < .Machine$double.eps] <- NA
        pc <- abs((mon[4:nrow(mon), (ncol(mon) - 1)] - prgam)/prgam)
        ifelse(pc[!is.na(pc)] > 0.01, TRUE, FALSE)
    }
    att <- 1
    while (!object$converge) {
        object <- update(object, step = 0.01, ...)
        att <- att + 1
        if (att > attempts) {
            mp.error.code("unstable")
            break
        }
    }
    att <- 1
    while (any(fwarn(object))) {
        object <- update(object, step = 0.01, ...)
        att <- att + 1
        if (att > attempts) {
            mp.error.code("converge")
            break
        }
    }
    object
}
