## Function pchisq.mixture
##
##' Probability distribution function for specified quantiles of a
##' mixture of chi-square distributions
##'
##' Probabilities for given quantiles are found from a mixture of
##' chi-square distribution based on mixing probabilities from a
##' binomial distribution with size equal to the number of traits or
##' environments and probability of "success" equal to 0.5.
##' @title Probabilities for specified quantiles of a mixture of
##' chi-square distributions
##' @param x Quantiles of the mixture distribution.
##' @param ntrait Number of traits or environments (or the size in
##' the underlying binomial distribution)
##' @return Probabilities for the specified quantiles are returned.
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##' @examples
##' ## Simple example
##' pchisq.mixture(5, ntrait=4)
##'
pchisq.mixture <- function(x, ntrait=2) {
  df <- 0:ntrait
  mixprobs <- dbinom(df, size=ntrait, prob=0.5)
  p <- c()
  for (i in 1:length(x)){
    p[i] <- sum(mixprobs*pchisq(x[i], df))
  }
  p
}
