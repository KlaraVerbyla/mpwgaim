## Function mpch.modify
##
## Internal function for modifying bivariate corh structures used in
## the models
##
## Calls \code{mpch.init} for terms in the model that involve
## bivariate corh structures
## @title Modification of initial estimates of corh
## structures
## @param model the model containing the corh structure
## @param Trait The trait or environment character vector
## @param by The genotype factor (present in the polygenic component
## of the model if it exists)
## @return The model structure is returned with modified initial
## estimates for the parameters of the corh components of
## the model
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mpch.modify <- function(model, Trait, by) {
  which.term <- igrep(list(Trait, 'ints'), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- mpch.init(Trait, 'ints', model$G.param)
  }
  which.term <- igrep(list(Trait, by), names(model$G.param))
  if(length(which.term) > 0) {
    model$G.param <- mpch.init(Trait, by, model$G.param)
  }
  model
}
