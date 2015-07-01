##  Functions for three-point probabilities
##
## These functions are not intended for direct use.  They are used in
## \code{mpcross2int}.
##
## xprobs(x, r13=0.1)
## xmult(x, r13=0.1)
## amult(x, r13=0.1)
## aprobs(x, r13=0.1)
## bprobs(x, r13=0.1)
## r2d(r=0.1)
## d2r(d=0.1)
## int4(x, which.prob=1, r13=0.1)
## int8(x, which.prob=1, r13=0.1)
## @title mpcross2int internal functions
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
xmult <- function(x, r13=0.1) {
    r23 <- (r13-x)/(1-2*x)
    xm <- c()
    xm[1] <- (1-x)*(1-r23)/2
    xm[2] <- (1-x)*r23/2
    xm[3] <- x*(1-r23)/2
    xm[4] <- x*r23/2
    xm[5] <- (1-x)/4
    xm[6] <- (1-r23)/4
    xm[7] <- (1-r13)/4
    xm[8] <- x/4
    xm[9] <- r23/4
    xm[10] <- r13/4
    xm
}
