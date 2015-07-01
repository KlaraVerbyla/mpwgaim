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
bprobs <- function(x, r13=0.1) {
    apattern <- c(rep(1,4),2:7,2:7,8:10)
    a.prob <- aprobs(x, r13)
    a.prob <- a.prob[apattern]
    a.mult <- amult(x, r13)
    b.prob <- a.prob*a.mult
    b.prob
}
