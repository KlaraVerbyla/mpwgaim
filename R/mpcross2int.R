## Function mpcross2int
##
##' Converts an object of class \code{"mpcross"} to an object with
##' class \code{"mpInterval"}.
##'
##' This function provides the conversion of genetic data object of
##' class \code{"mpcross"} from the \pkg{mpMap} package to an
##' \code{"mpInterval"} object that is ready to use with \pkg{mpwgaim}
##' QTL analysis functions.
##'
##' @title Convert an mpcross object to an mpInterval object
##' @param mpcross An \code{mpcross} object generated for example
##' using the \pkg{mpMap} package.
##' @param gen.type Founder probabilities found for intervals on the
##' linkage map using \code{"mpInterval"} or markers using
##' \code{"mpMarker"}.  Default is \code{"mpMarker"}.
##' @param method Method for calculation of founder probabilities.
##' Default is \code{"qtl"} to use Hidden Markov Models in the
##' \pkg{qtl} package.  The alternative is \code{"3pt"} to use
##' probabilities found using haplotypes of length 3.
##' @param step Step size for calculation of founder probabilities
##' based on the Hidden Markov Model using integration (Trapezoidal
##' rule)
##' @param population A character vector which indicates the numer of
##' founders in the MAGIC population, required for the calculation of
##' three point probabilities.  The default is \code{"4way"} (the
##' other current option is \code{"8way"}).
##' @param debug Logical (default is \code{FALSE}) to determine if
##' printing is to be enabled in debugging the function
##' @return An mpInterval object (list) suitable for QTL analysis using
##' mpwgaim.  The list contains the following components:
##' \itemize{
##' \item \code{geno} : A list with elements named after the
##' corresponding names of the chromosomes or linkage groups.  Each
##' component of the list has potentially five components.  These are
##' the linkage map in \code{map}, the marker data \code{data}, marker
##' data for the founders in \code{founders}, and \code{dist} which
##' has the distances in cM for the linkage group.  The fifth
##' component depends on whether the analysis using \code{mpwgaim} is
##' based on intervals (component \code{intval} or markers (component
##' \code{imputed.data}).
##' \item \code{pheno} : A vector of four-way line names
##' \item \code{founders} : A vector of founder names
##' \item \code{nfounders} : The number of founders in the MAGIC
##' population
##' }
##' @author Ari Verbyla (ari.verbyla at csiro.au)
##' @export
##' @examples
##' \dontrun{
##' ## Simulated data to illustrate generating the interval object.
##' ## We require the mpMap and qtl packages:
##' require(mpMap)
##' require(qtl)
##' set.seed(3)
##' map <- sim.map(len=rep(300,7), n.mar=rep(201,7), eq.spacing=TRUE,
##'                include.x=FALSE)
##' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
##' ## Dummy QTL to set up structure
##' qtl.mat <- matrix(data=c(1, 142, 0.3, -0.3, -0.3, 0.3,
##'                          2, 162, 0.3, -0.3, -0.3, 0.3,
##'                          3, 78, 0, 0, 0, 0,
##'                          4, 114, 0.3, -0.3, -0.3, 0.3),
##'                   nrow=4, ncol=6, byrow=TRUE)
##' ## Set up the mpcross object
##' sim.dat <- sim.mpcross(map=map,
##' pedigree=sim.ped, qtl=qtl.mat, seed=5)
##' mpSim <- maporder(sim.dat)
##' ## Population is four way which is the default setting.
##' ## Interval object for analysis using intervals and Hidden Markov
##' ## model for probability calculations
##' mpInterval <- mpcross2int(mpSim, gen.type="mpInterval")
##' ## Interval object for analysis using markers and Hidden Markov
##' ## model for probability calculations
##' mpMarker <- mpcross2int(mpSim, gen.type="mpMarker")
##' ## Interval object for analysis using intervals and three-point probabilities
##' mpInterval <- mpcross2int(mpSim, gen.type="mpInterval", method="3pt")
##' ## Interval object for analysis using markers and three-point probabilities
##' mpMarker <-mpcross2int(mpSim, gen.type="mpMarker", method="3pt")
##' }
##'
mpcross2int <- function(mpcross, gen.type="mpMarker", method = "qtl",
                        step = 0.1, population = "4way", debug=FALSE) {
    if(!any(class(mpcross) == "mpcross"))
        stop(" Argument is not an mpcross object\n")
    if(!gen.type %in% c("mpMarker", "mpInterval"))
        stop("gen.type needs to be either mpMarker or mpInterval")
    nfounders <- dim(mpcross$founders)[1]
    ngeno <- length(mpcross$map)
##    require(mpMap)  delete - ABZ
    ##    mpcross <- maporder(mpcross)
    ##    mpcross <- newcleanmap(mpcross, mindist=0.1)
    ##
    magicObj <- list()
    magicObj$pheno <- data.frame(id=rownames(mpcross$finals))
    magicObj$geno <-  vector("list", ngeno)
    magicObj$nfounders <- nfounders
    magicObj$founders <- dimnames(mpcross$founders)[[1]]
    names.geno <- names(magicObj$geno) <- names(mpcross$map)
    nmark <- unlist(lapply(mpcross$map, function(x) length(x)))
    markers <- mpcross$finals
    founders <- mpcross$founders
    n.alleles <- apply(markers, 2, function(x) length(unique(x[!is.na(x)])))
    nonbi <- n.alleles > 2
    tmarkers <- matrix(as.character(markers), ncol=ncol(markers))
    tfounders <- matrix(as.character(founders), ncol=ncol(founders))
    for (ii in 1:dim(tmarkers)[2]) {
        if(!nonbi[ii]) {
            uni.ii <- sort(unique(tmarkers[,ii]))
            uni.ii <- uni.ii[!is.na(uni.ii)]
            tmarkers[, ii] <- gsub(uni.ii[1], -1, tmarkers[,ii])
            tmarkers[, ii] <- gsub(uni.ii[2], 1, tmarkers[,ii])
            tfounders[, ii] <- gsub(uni.ii[1], -1, tfounders[,ii])
            tfounders[, ii] <- gsub(uni.ii[2], 1, tfounders[,ii])
        }
    }
    new.markers <- matrix(as.numeric(tmarkers), ncol=ncol(markers))
    dimnames(new.markers) <- dimnames(markers)
    new.founders <- matrix(as.numeric(tfounders), ncol=ncol(founders))
    dimnames(new.founders) <- dimnames(founders)
    upper <- 0
    for (ii in 1:ngeno) {
        lower <- upper + 1
        upper <- upper + nmark[ii]
        magicObj$geno[[ii]]$map <- mpcross$map[[ii]]
        magicObj$geno[[ii]]$data <- as.matrix(new.markers[,lower:upper])
        magicObj$geno[[ii]]$founders <- as.matrix(new.founders[,lower:upper])
        dimnames(magicObj$geno[[ii]]$founders)[[2]] <- paste("Chr", names(magicObj$geno)[ii],
                                                             1:nmark[ii], sep='.')
        magicObj$geno[[ii]]$dist <- magicObj$geno[[ii]]$map
    }
    ##
    tmp <- vector(mode="list", length=ngeno)
    if(gen.type == "mpInterval") {
        if(method == "3pt") {
            cat(" Creating interval object using 3-point haplotype probabilities\n")
            if(population == "4way") lookup <- lookup4
            else lookup <- lookup8
            intfunc <- paste("int", nfounders, sep="")
            nprobs <- 10
            if(nfounders == 8) nprobs <- 19
            for (cc in 1:length(magicObj$geno)) {
                cat(" linkage group ", cc, "\n")
                dset <- diff(magicObj$geno[[cc]]$map)/100
                lg.prob <- c()
                for (ii in 1:length(dset)) {
                    if(debug) print(ii)
                    r0 <- d2r(d=dset[ii])
                    int.prob <- matrix(0,1,nfounders)
                    ab.prob <- c()
                    for (aa in 1:nprobs) {
                        ab.prob[aa] <- integrate(Vectorize(eval(parse(text=intfunc))), 0, r0,
                                                 which.prob = aa, r13=r0)$value
                    }
                    if(debug) print(ab.prob)
                    for (ll in 1:dim(magicObj$geno[[cc]]$data)[1]) {
                        ML <- magicObj$geno[[cc]]$data[ll,ii]
                        if(debug) print(ML)
                        MR <- magicObj$geno[[cc]]$data[ll,ii+1]
                        if(debug) print(MR)
                        if(is.na(ML)) left.set <- 1:nfounders
                        else left.set <- (1:nfounders)[magicObj$geno[[cc]]$founders[,ii] == ML]
                        if(debug) {
                            print(magicObj$geno[[cc]]$founders[,ii])
                            print(left.set)
                        }
                        if(is.na(MR)) right.set <- 1:nfounders
                        else right.set <- (1:nfounders)[magicObj$geno[[cc]]$founders[,ii+1] == MR]
                        if(debug) print(right.set)
                        Qprob <- rep(1/nfounders, nfounders)
                        if(!(is.na(ML) & is.na(MR))) {
                            Qprob <- rep(0,nfounders)
                            for (ls in left.set) {
                                for (rs in right.set) {
                                    for (qq in 1:nfounders) {
                                        hap.pattern <- paste(ls,qq,rs,sep='')
                                        if(debug) cat(" pattern ", hap.pattern, "\n")
                                        which.prob <- lookup$pattern[lookup$haplotype == hap.pattern]
                                        if(debug) cat(" which.prob ", which.prob, " ab.prob ", ab.prob[which.prob], "\n")
                                        Qprob[qq] <- Qprob[qq] + ab.prob[which.prob]
                                    }
                                }
                            }
                            Qprob <- Qprob/sum(Qprob)
                        }
                        int.prob <- rbind(int.prob, Qprob)
                    }
                    int.prob <- int.prob[-1,]
                    lg.prob <- cbind(lg.prob, int.prob)
                }
                tmp[[cc]] <- lg.prob
                dimnames(tmp[[cc]])[[1]] <- paste("L", 1:dim(magicObj$geno[[cc]]$data)[1],sep="")
            }
        }
        else {
            cat(" Creating interval object using HMM haplotype probabilities\n")
            cc.count <- 0
            for (cc in names(mpcross$map)) {
                cc.count <- cc.count+1
                cat(" Linkage group ", cc, "\n")
                exp.prob <- mpprob(mpcross, chr = cc, program="qtl", step=step)$prob
                if(!class(dimnames(exp.prob[[cc]])[[2]]) == "character") dimnames(exp.prob[[cc]])[[2]] <- as.character(dimnames(exp.prob[[cc]])[[2]])
                name.prob <- unlist(lapply(strsplit(dimnames(exp.prob[[cc]])[[2]], split=","), function(el) el[1]))
                spl.prob <- unique(name.prob)
                which.loc <- substr(spl.prob, start=1, stop=3) == "loc"
                loc <- vector("numeric", length=length(spl.prob))
                spl.prob1 <- spl.prob[which.loc]
                loc1 <- unlist(lapply(strsplit(spl.prob1, split="c"), function(el) as.numeric(el[2])))
                loc[which.loc] <- loc1
                loc[!which.loc] <- mpcross$map[[cc.count]]
                new.probs <- c()
                for (ii in 1:(length(mpcross$map[[cc.count]])-1) ) {
                    mname1 <- names(mpcross$map[[cc.count]])[ii]
                    if(debug) print(mname1)
                    mname2 <- names(mpcross$map[[cc.count]])[ii+1]
                    if(debug) print(mname2)
                    which.p1 <- grep(mname1, name.prob, fixed=TRUE)
                    if(debug) print(which.p1)
                    which.p1 <- which.p1[name.prob[which.p1] %in% mname1]
                    if(debug) print(which.p1)
                    which.p2 <- grep(mname2, name.prob, fixed=TRUE)
                    if(debug) print(which.p2)
                    which.p2 <- which.p2[name.prob[which.p2] %in% mname2]
                    if(debug) print(which.p2)
                    int.probs <- exp.prob[[cc]][,min(which.p1):max(which.p2)]
                    if(debug) print(dim(int.probs))
                    which.m1 <- grep(mname1, spl.prob, fixed=TRUE)
                    if(debug) print(which.m1)
                    which.m1 <- which.m1[spl.prob[which.m1] %in% mname1]
                    if(debug) print(which.m1)
                    which.m2 <- grep(mname2, spl.prob, fixed=TRUE)
                    if(debug) print(which.m2)
                    which.m2 <- which.m2[spl.prob[which.m2] %in% mname2]
                    if(debug) print(which.m2)
                    dseq <- loc[which.m1:which.m2]
                    if(debug) print(dseq)
                    dseq <- diff(dseq)
                    if(debug) print(dseq)
                    for (ff in 1:nfounders) {
                        if(debug) print(ff)
                        seqf <- seq(from=which.p1[ff], to=which.p2[ff], by=nfounders)-min(which.p1)+1
                        probs1 <- int.probs[, seqf]
                        if(debug) print(dim(probs1))
                        ff.probs <- apply(probs1, 1, function(el, dseq)
                                          sum((el[-length(el)]+el[-1])*dseq)/(2*(mpcross$map[[cc.count]][ii+1]-
                                                                                 mpcross$map[[cc.count]][ii])), dseq)
                        if(any(ff.probs>1 | ff.probs < 0)) {
                            cat(" Probabilities outside [0,1] on linkage group ", cc, " interval ", ii,
                                " founder ", ff, "\n")
                            stop()
                        }
                        if(any(is.na(ff.probs))) {
                            cat("  Missing probabilities, linkage group ", cc, " interval ", ii,
                                " founder ", ff, "\n")
                            if(debug) print(sum(is.na(ff.probs)))
                        }
                        new.probs <- cbind(new.probs, ff.probs)
                    }
                }
                tmp[[cc.count]] <- new.probs
            }
        }
    }
    else {
        if(method == "3pt") method <- "mpMap"
        tmp <- mpprob(mpcross, program=method, step=0)$prob
    }
    ##
    for (ii in 1:ngeno) {
        chrnam <- names.geno[ii]
        lchr <- dim(tmp[[ii]])[2]/nfounders
        tmp1 <- paste('Chr', chrnam, 1:lchr, sep='.')
        dimnames(tmp[[ii]])[[2]] <- paste(rep(tmp1,each=nfounders),
                                          rep(1:nfounders,lchr), sep='.')
    }
    if(gen.type == "mpInterval") {
        for (ii in 1:ngeno) {
            magicObj$geno[[ii]]$intval <- tmp[[ii]]
        }
    }
    else if(gen.type == "mpMarker") {
        for (ii in 1:ngeno) {
            magicObj$geno[[ii]]$imputed.data <- tmp[[ii]]
        }
    }
    ##
    magicObj$gen.type <- gen.type
    class(magicObj) <- c("mpInterval", "cross", "interval")
    magicObj
}
