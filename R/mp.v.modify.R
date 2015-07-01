##  Internal starting value functions
##
## Internal functions for resetting initial estimates of variance
## components
##
## mp.v.modify(model)
## mp.v.init(term, G.param)
## @title Internal functions to modify zero variance components
## @param model model with marker/interval effects and QTL effects
## @return These functions reset zero variance components to a small
## value to ensure 'asreml' fits the components in a non-fixed state.
## @author Ari Verbyla (ari.verbyla at csiro.au)
##
mp.v.modify <- function(model)
{
    which.term <- names(model$G.param)[grep("X\\.", names(model$G.param))]
    if (length(which.term) > 0) {
        model$G.param <- mp.v.init(which.term, model$G.param)
    }
    which.term <- grep("grp(\"ints\")", names(model$G.param))
    if (length(which.term) > 0) {
        model$G.param <- mp.v.init("grp(\"ints\")", model$G.param)
    }
    model
}
