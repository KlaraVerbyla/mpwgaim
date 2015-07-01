## Function mp.envDestruct
##
## Internal function to minimise the size of 'mpwgaim' objects
##
## Used in 'mpwgaim.asreml' on exit to reduce the size of the final
## object.
## @title Reduction of information in an 'mpwgaim' fit.
## @param model Final model in 'mpwgaim'.
## @param keep Character vector of objects to keep.
## @return A reduced model object.
## @author Julian Taylor (julian.taylor at adelaide.edu.au)
##
mp.envDestruct <- function(model, keep)
{
    fixed <- deparse(model$fixed.formula)
    random <- deparse(model$random.formula)
    if (!is.null(sparse <- model$call$sparse)) {
        if ("mv" %in% all.vars(sparse))
            sparse <- update.formula(sparse, ~. - mv)
        if (is.numeric(sparse[[2]]))
            sparse <- NULL
        else sparse <- deparse(sparse)
    }
    if ("rcov" %in% names(model$call))
        rcov <- deparse(model$call$rcov)
    lsp <- ls(parent.frame(n = 1))
    rm(list = lsp[!(lsp %in% keep)], pos = parent.frame(n = 1))
    model$call$fixed <- formula(fixed)
    model$call$random <- formula(random)
    model$fixed.formula <- parse(text = fixed)[[1]]
    model$random.formula <- parse(text = random)[[1]]
    if (!is.null(sparse)) {
        model$call$sparse <- formula(sparse)
        model$sparse.formula <- parse(text = sparse)[[1]]
    }
    if ("rcov" %in% names(model$call))
        rcov <- formula(rcov)
    model
}
