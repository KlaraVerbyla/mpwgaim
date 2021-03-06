% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/mpwgaim.R
\name{mpwgaim}
\alias{mpwgaim}
\title{mpwgaim method for class \code{"asreml"}}
\usage{
mpwgaim(baseModel, ...)
}
\arguments{
\item{baseModel}{a model object of class "'asreml'" usually
representing a base model with which to build the qtl model.}

\item{\ldots}{Further arguments to be passed to method
\code{\link{mpwgaim.asreml}}.}
}
\value{
An object of class "'mpwgaim'" which also inherits the
class "'asreml'" by default. The object returned is actually an
'asreml' object (see 'asreml.object') with the addition of
components from the QTL detection listed below.  \item{QTL:}{ A
list of components from the significant QTL detected including a
character vector of the significant QTL along with a vector of the
QTL effect sizes. There are also a number of diagnostics that can
be found in the 'diag' component.}\
}
\description{
Fits an iterative Multi-Parent Whole Genome Average Interval
Mapping (mpwgaim) model for QTL detection
}
\details{
This function implements the whole genome average interval mapping
approach for multi-parent populations.  It includes code to handle
high-dimensional problems.  The size of computations in model
fitting is the number of genetic lines rather than the number of
markers or intervals.  The use of mpwgaim requires a base model to
be fitted which includes all genetic and non-genetic terms that
are relevant for an analysis without marker information.  The
marker information must be contained in an interval object
('intervalObj') and these effects are added to the phenotypic data
('phenoData') to allow fitting of models with QTL.  Note that
analysis using 'markers' or 'intervals' can be carried out
('gen.type' is set by default to 'interval').  All QTL selected
are random effects in the final model.
}
\examples{
## Simulated data example
\dontrun{
require(asreml)
require(mpMap)
require(qtl)
map <- sim.map(len=rep(200,7), n.mar=rep(51,7), eq.spacing=TRUE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
## QTL to set up mpcross structure
qtl.mat <- matrix(data=c(1, 142, 0.354, -0.354, -0.354, 0.354,
                  2, 162, 0.354, -0.354, -0.354, -0.354,
                  5, 78, 0.354, -0.354, -0.354, 0.354),
                  nrow=3, ncol=6, byrow=TRUE)
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped,
                       qtl=qtl.mat, seed=5)
mpSim <- maporder(sim.dat)
## Use Hidden Markov model for probability calculations
mpInterval <- mpcross2int(mpSim, gen.type="mpInterval")
mpMarker <- mpcross2int(mpSim, gen.type="mpMarker")
##  Model matrix for computing the contribution of the QTL
nqtl.col <- dim(sim.dat$qtlgeno$finals)[2]
mmat <- matrix(nrow=500, ncol=4*nqtl.col/2)
for (ii in 1:(nqtl.col/2)) {
  qtl.fac <- factor(sim.dat$qtlgeno$finals[,ii])
  mmat[,(4*ii-3):(4*ii)] <- model.matrix(~qtl.fac - 1)
}
## Effects
qtl.sizes <- as.vector(t(qtl.mat[,3:6]))
qtl.effect <- mmat \%*\% qtl.sizes
## Polygenic variance
pvar <- 0.5
## Function to calculate approximate percentage variance for each QTL
perc.var <- function(qtl.mat, poly.var) {
    nfounders <- dim(qtl.mat)[2]-2
    prob <- 1/nfounders
    varq <- diag(rep(prob,nfounders)) - rep(prob,nfounders) \%*\% t(rep(prob,nfounders))
    gvar <- apply(qtl.mat[, -c(1,2)], 1, function(el, varq) sum((el \%*\% varq) * el), varq)
    totvar <- sum(gvar)+pvar
    perc.var <- 100*gvar/totvar
    round(perc.var,1)
}
## Percentage variance for each QTL
percvar <- perc.var(qtl.mat, pvar)
percvar
## Setup simulated data for analysis
ngeno <- 500
nrep <- 2
id <- factor(rep(paste0("L", 1:ngeno),nrep))
pheno.data <- data.frame(y = 0, Rep=nrep, id=id)
set.seed(10)
ee <- rnorm(ngeno*nrep,0,1)
uu <- rnorm(ngeno,0,1)
pheno.data$y <- 10 + rep(qtl.effect,2) + sqrt(1/2)*c(uu,uu) + ee
sim.asr0 <- asreml(y ~ 1, random = ~ id, data=pheno.data)
sim.qtl <- mpwgaim(sim.asr0, pheno.data, mpInterval, merge.by = "id",
                  verboseLev=0, gen.type="interval", na.method.X='include',
                  data.name = "sim.data")
sim.summ <- summary(sim.qtl, mpInterval)
sim.summ
plot(sim.summ, mpInterval)
## Can repeat the analysis using mpMarker instead of mpInterval
}
}
\author{
Ari Verbyla (ari.verbyla at csiro.au) and Klara Verbyla
(klara.verbyla at csiro.au)
}
\references{
Verbyla, A. P., George, A. W., Cavanagh, C. C. and
Verbyla, K. L. (2014).  Whole genome QTL analysis for MAGIC.
Theoretical and Applied Genetics.  DOI:10.1007/s00122-014-2337-4.
}
\seealso{
'print.mpwgaim', 'summary.mpwgaim, 'plot.summary.mpwgaim'
}

