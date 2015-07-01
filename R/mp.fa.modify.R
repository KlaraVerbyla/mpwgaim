## Function mp.fa.modify
##
## Internal function for modifying factor analytic structures used in
## the models
##
## Calls \code{mp.favar.init} for terms in the model that involve
## factor analytic structures
## @title Modification of initial estimates of factor analytic
## structures
## @param model the model containing the factor analytic structure
## @param Trait The trait or environment character vector
## @param by The genotype factor (present in the polygenic component
## of the model if it exists)
## @return The model structure is returned with modified initial
## estimates for the parameters of the factor analytic components of
## the model
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mp.fa.modify <- function(model, Trait, by) {
  which.term <- igrep(list(paste('fa\\(', Trait, sep=''), 'ints'), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- mp.favar.init(Trait, 'ints', model$G.param)
  }
  which.term <- igrep(list(paste('fa\\(', Trait, sep=''), by), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- mp.favar.init(Trait, by, model$G.param)
  }
  model
}
