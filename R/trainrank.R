#' R wrapper for the train_rank.py script.
#'
#' For larger numbers of metabolites (>50) calling this function will
#' take a while.
#' @title The metabolites/transcripts HMM association program
#' @param mdat matrix with the metabolite data (time points in columns)
#' @param tdat matrix with the transcript data (time points in columns)
#' @param states number of states
#' @param cyclic use cyclic model
#' @param LD_LIBRARY_PATH path to ghmm shared object files
#' @param trainrank where to find the train_rank.py script
#' @param python the command use to call python
#' @param verbose be verbose
#' @param cleanOnExit if true, remove temp files and created
#' directories on exit
#' @param scale logical indicating if the data should be scaled to
#' unitvariance and centered. HMMs are very sensitive to this and
#' scaling is strongly recommended.
#' @return a list with results for negative and positive modes
#' @references Redestig, H. and Costa, IG (submitted)
#' @examples
#' trainrank(rbind(sample(1:10)), rbind(rnorm(10)))
#' metabolite <- jitter(c(1,1,1,2,2,2,6,6,6,6,6,2,2,2))
#' associated_gene <- jitter(c(2,2,2,6,6,6,6,6,2,2,2,3,3,3))
#' unassociated_genes <- t(replicate(50, rnorm(14)))
#' mdat <- rbind(metabolite)
#' tdat <- rbind(associated_gene, unassociated_genes)
#' rownames(mdat) <- rownames(tdat) <- NULL
#' matplot(cbind(metabolite, associated_gene), type='b')
#' tt <- trainrank(mdat, tdat)
#' barplot(tt$positive[,1])
#' tt$plag[1,1]
#' @author Henning Redestig
#' @export
trainrank <- function(mdat, tdat, states=3, cyclic=FALSE,
                      LD_LIBRARY_PATH="/usr/local/lib/",
                      trainrank=system.file("extdata", "train_rank.py", package="trainmet"),
                      python="python -W ignore",
                      verbose=FALSE, cleanOnExit=TRUE, scale=TRUE) {
    hackm <- FALSE
    hackt <- FALSE
    if(nrow(mdat) < 2) {
        mdat <- rbind(mdat, mdat)
        rownames(mdat) <- make.unique(make.names(rownames(mdat)))
        hackm <- TRUE
    } 
    if(nrow(tdat) < 2) {
        tdat <- rbind(tdat, tdat)
        rownames(tdat) <- make.unique(make.names(rownames(tdat)))
        hackt <- TRUE
    }
    if(scale){
        mdat <- t(scale(t(mdat)))
        tdat <- t(scale(t(tdat)))
    }

    minfile <- tempfile()
    if(verbose) {
        cat("metabolites infile ")
        print(minfile)
    }
    tinfile <- tempfile()
    if(verbose) {
        cat("transcripts infile ")
        print(tinfile)
    }
    ## trainrank doesnt like output in other directory
    outfile <- tempfile(tmpdir=".")
    if(verbose) {
        cat("outfile ")
        print(outfile)
    }
    if(cleanOnExit) {
        on.exit(unlink(minfile))
        on.exit(unlink(tinfile), add=TRUE)
        on.exit(unlink(paste(outfile, ".res", sep="")), add=TRUE)
        on.exit(unlink(paste(minfile, ".sqd", sep="")), add=TRUE)
        on.exit(unlink(paste(tinfile, ".sqd", sep="")), add=TRUE)
        on.exit(unlink(paste(outfile, ".lag", sep="")), add=TRUE)
        on.exit(unlink(paste(outfile, "_n.lag", sep="")), add=TRUE)
        on.exit(unlink(paste(outfile, "_n.res", sep="")), add=TRUE)
    }
    if(!is.null(LD_LIBRARY_PATH))
        if(length(grep(LD_LIBRARY_PATH, Sys.getenv("LD_LIBRARY_PATH"))) == 0)
            Sys.setenv(LD_LIBRARY_PATH=paste(LD_LIBRARY_PATH,
                           Sys.getenv("LD_LIBRARY_PATH"), sep=":"))

    write.table(cbind(1:nrow(mdat),mdat),
                file=minfile, row.names=TRUE, quote=FALSE, col.names=FALSE,
                sep="\t")
    write.table(cbind(1:nrow(tdat),tdat),
                file=tinfile, row.names=TRUE, quote=FALSE, col.names=FALSE,
                sep="\t")

    if(!file.exists("models")) {
        dir.create("models")
        if(cleanOnExit) on.exit(unlink("models", recursive=TRUE), add=TRUE)
    }

    system(paste(python, trainrank, ifelse(cyclic, "-c", ""),
                 "-s", states, minfile, tinfile, outfile,
                 ifelse(verbose, "", "> /dev/null")))
    pos <- read.table(paste(outfile, ".res", sep=""), sep="\t", quote="",
                      header=TRUE, row.names=1, check.names=FALSE)
    neg <- read.table(paste(outfile, "_n.res", sep=""), sep="\t", quote="",
                      header=TRUE, row.names=1, check.names=FALSE)  
    poslag <- read.table(paste(outfile, ".lag", sep=""), sep="\t", quote="",
                         header=TRUE, row.names=1, check.names=FALSE)
    neglag <- read.table(paste(outfile, "_n.lag", sep=""), sep="\t", quote="",
                         header=TRUE, row.names=1, check.names=FALSE)  

    if(hackm) {
        pos <- pos[,1,drop=FALSE]
        neg <- neg[,1,drop=FALSE]
        poslag <- poslag[,1,drop=FALSE]
        neglag <- neglag[,1,drop=FALSE]

    }
    if(hackt) {
        pos <- pos[1,,drop=FALSE]
        neg <- neg[1,,drop=FALSE]
        poslag <- poslag[1,,drop=FALSE]
        neglag <- neglag[1,,drop=FALSE]
    }    
    list(positive=pos, negative=neg, plag=poslag, nlag=neglag)
}

