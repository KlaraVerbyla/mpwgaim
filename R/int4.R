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
int4 <- function(x, which.prob=1, r13=0.1) {
    a.prob <- aprobs(x, r13=r13)
    a.prob <- a.prob[which.prob]
    dist <- r2d(r=r13)
    divisor <- dist*(1-2*x)
    aint <- a.prob/divisor
    aint
}
