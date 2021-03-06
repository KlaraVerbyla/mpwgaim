% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/summary.mpwgaim.R
\name{summary.mpwgaim}
\alias{summary.mpwgaim}
\title{Summary method for the class 'mpwgaim'}
\usage{
\method{summary}{mpwgaim}(object, intervalObj, merge.by = "id", ...)
}
\arguments{
\item{object}{an object of class 'mpwgaim'}

\item{intervalObj}{an object of class 'interval'}

\item{merge.by}{The genetic line variable used in merging
phenotypic and genetic data in the 'mpwgaim' fit.}

\item{\ldots}{Optional additional objects passed to or from other
methods}
}
\value{
A object of class 'summary.mpwgaim, consisting of a list
with a data frame of QTL information to be used with the print
method, components 'LODtrace', 'qtl' and 'type' that are used for
the plotting method 'plot.summary.mpwgaim'
}
\description{
Prints a summary of the "'mpwgaim'" object in a presentable format
}
\details{
After a QTL analysis using 'mpwgaim', a summary of the QTL can be
obtained using this function.  In addition the building blocks are
calculated for a plot of log-probabilities across the genome
including QTL.
}
\examples{
\dontrun{
## Simulated data example
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
\seealso{
'print.summary.mpwgaim', 'plot.summary.mpwgaim'
}

