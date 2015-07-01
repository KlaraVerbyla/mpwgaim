## Function qchisq.mixture
##
##' Internal function for quantiles of a mixture of chi-square
##' distributions
##'
##' An iterative process is used to determine the quantiles for
##' specified probabilities in \code{prob}, with mixing probabilities in
##' the mixture from a binomial distribution with size the number of
##' traits or environments (in the context of mvmpwgaim) and
##' probability of success equal to 0.5.  The process is stable for
##' probabilities that are not too small.  For a critical value in a
##' test using this distribution with type I error 0.05, the
##' probability should be specified as 0.95 (so that the probability is
##' assumed to relate to values of the mixture distribution less than
##' the critical value).
##' @title Quantiles of a mixture of chi-square distributions
##' @param prob A single or vector of probabilities (
##' @param ntrait Number of traits or environments (or the binomial
##' size)
##' @param trace Logical set to \code{FALSE}.  If \code{TRUE} then
##' the iterative process is printed.
##' @param maxiter Maximum number of iterations to carry out.  The
##' default is 10.
##' @return The quantiles for the specified probabilities are returned.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##' @examples
##' ## Simple examples
##' qchisq.mixture(0.95, ntrait=4)
##' qchisq.mixture(0.05, ntrait=4) ## fails
##'
qchisq.mixture <- function(prob, ntrait=2, trace=F, maxiter=10) {
  ##
  #   Iterative procedure to calculate percentiles for mixtures of chisquare variates
  #
  #   df the degrees of freedom for the constituent chi-squared variates
  #   mixprobs  the mixing probabilities
  #
  df <- 0:ntrait
  mixprobs <- dbinom(df, size=ntrait, prob=0.5)
  ##
  #  Start value
  ##
  cv <- qchisq(prob, df=ntrait)
  if(trace) cat(" Starting value: ", cv, "\n")
  obj.fn <- function(cv=cv, df=cv, mixprobs=mixprobs) {
    obj <- ifelse(df[1] == 0, mixprobs[1] + sum(pchisq(cv, df=df[-1])*mixprobs[-1]),
                  sum(pchisq(cv, df=df)*mixprobs)) - prob
    d.obj <- ifelse(df[1] == 0, sum(dchisq(cv, df=df[-1])*mixprobs[-1]),
                    sum(dchisq(cv, df=df)*mixprobs))
    list(f=obj, df=d.obj)
  }
  convergence <- F
  iteration <- 0
  while(!convergence & iteration <= maxiter) {
    obj <- obj.fn(cv=cv, df=df, mixprobs=mixprobs)
    corr <- (obj$f/obj$df)
    cv <- cv - corr
    iteration <- iteration + 1
    if(trace) cat(" Iteration ", iteration, ": cv = ", cv, "\n")
    convergence <- abs(corr) < 0.00001
  }
  if(!convergence) warning(' non-convergence: result incorrect\n')
  cv
}
