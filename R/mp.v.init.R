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
mp.v.init <- function(term, G.param)
{
  for (i in 1:length(term)) {
        which.term <- grep(term[i], names(G.param), fixed=TRUE)
        which.term <- which.term[names(G.param)[which.term] %in%
            term[i]]
        var.term <- grep("var", names(G.param[[which.term]][[1]]$initial))
        con.term <- G.param[[which.term]][[1]]$con[var.term] ==
            "B"
        if (con.term) {
            G.param[[which.term]][[1]]$con[con.term] <- "P"
            G.param[[which.term]][[1]]$initial[con.term] <- 0.1
        }
    }
    G.param
}
