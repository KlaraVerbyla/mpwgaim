## Function mpch.init
##
## Internal function to initialise the specific variance parameters
## that have been fixed in fitting corh structures
##
## Those components of a fitted corh model that have been
## fixed in the fitting process are reset for subsequent fitting.
## @title Modification of initial estimates for fixed components of
## corh structures
## @param char1 The name of the trait or environment factor as a
## character
## @param char2 A character vector, either '"ints'" or the genotype
## name
## @param G.param the G.param list of a fitted asreml model
## @return A new asreml G.param list that changes the correlation
## that has been fixed close to 1 in the fitting of the corh
## model.
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mpch.init <- function(char1, char2, G.param) {
  which.term <- igrep(list(char1, char2), names(G.param))
  which.subterm <- grep(char1, names(G.param[[which.term]]))
  cor.term <- grep('cor', names(G.param[[which.term]][[which.subterm]]$initial))
#  print(cor.term)
  con.term <- G.param[[which.term]][[which.subterm]]$con[cor.term] == 'B'
#  print(con.term)
  if(sum(con.term) > 0) {
    G.param[[which.term]][[which.subterm]]$con[cor.term] <- 'U'
    G.param[[which.term]][[which.subterm]]$initial[cor.term] <- 0.1
  }
  G.param
}
