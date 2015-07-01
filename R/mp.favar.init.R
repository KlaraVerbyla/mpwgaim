## Function mp.favar.init
##
## Internal function to initialise the specific variance parameters
## that have been fixed in fitting factor analytic structures
##
## Those components of a fitted factor analytic model that have been
## fixed in the fitting process are reset for subsequent fitting.
## @title Modification of initial estimates for fixed components of
## factor analytic structures
## @param char1 The name of the trait or environment factor as a
## character
## @param char2 A character vector, either '"ints'" or the genotype
## name
## @param G.param the G.param list of a fitted asreml model
## @return A new asreml G.param list that changes those specific
## variances that have been fixed at zero in the fitting of the
## model.
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mp.favar.init <- function(char1, char2, G.param) {
  which.term <- igrep(list(char1, char2), names(G.param))
  which.subterm <- grep(char1, names(G.param[[which.term]]))
  var.terms <- grep('var', names(G.param[[which.term]][[which.subterm]]$initial))
  con.terms <- G.param[[which.term]][[which.subterm]]$con[var.terms] == 'F'
  if(sum(con.terms) > 0) {
    G.param[[which.term]][[which.subterm]]$con[con.terms] <- 'P'
    G.param[[which.term]][[which.subterm]]$initial[con.terms] <- 0.001
  }
  G.param
}
