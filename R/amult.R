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
amult <- function(x, r13=0.1) {
    r23 <- (r13-x)/(1-2*x)
    am <- c()
    am[1] <- (1-x)*(1-r23)/2
    am[2] <- (1-x)*r23/2
    am[3] <- x*(1-r23)/2
    am[4] <- x*r23/2
    am[5] <- (1-x)/4
    am[6] <- (1-r23)/4
    am[7] <- (1-r13)/4
    am[8] <- (1-x)/4
    am[9] <- (1-r23)/4
    am[10] <- (1-r13)/4
    am[11] <- x/4
    am[12] <- r23/4
    am[13] <- r13/4
    am[14] <- x/4
    am[15] <- r23/4
    am[16] <- r13/4
    am[17] <- am[18] <- am[19] <- 1/8
    am
}
