##' Labels are made by taking the target columns from fData and
##' filling labels sequentially from the columns in fData prioritizing
##' the selected columns. Labels are then made unique and any
##' remaining missing values are replaced by a constructed identifier.
##' @title Make unique non-mising labels for a set of metabolites
##' @param object an ExpressionSet
##' @param targets the selected columns in fData to search in.
##' @return character vector with non-missing unique elements
##' @author henning
##' @examples
##' data(mdat)
##' makeMetLabels2(mdat)
##' @export
makeMetLabels2 <- function(object, targets=c("preferred","synonym", "lineno")) {
  if(!any(targets %in% names(fData(object))))
    stop("none of the requested columns in the fData")
  targets <- c(targets, names(fData(object))[!names(fData(object)) %in% targets])
  res <- fData(object)[,targets[1]]
  res[res == ""] <- NA
  for(i in 2:length(targets)) {
    if(!any(is.na(res)))
      break
    res[is.na(res)] <- fData(object)[,targets[i]][is.na(res)]
    res[res == ""] <- NA
  }
  if(any(is.na(res)))
    res[is.na(res)] <- paste("NA_", 1:sum(is.na(res)), sep="")
  make.unique(res)
}

##' Make a line plot of transcription and metabolite abundnance
##' trajectories over a time-course and color the genes after their
##' ranked co-response to the metabolite.
##'
##' Can be used to visualize how highly each metabolite ranks the
##' different transcripts. Ranks are relative to the all genes but
##' only the selected genes are plotted.
##' @title Color line plot
##' @param tdat the gene expression data. Typically the whole dataset.
##' @param mdat the metabolite data. Only the metabolites of
##' interest. 
##' @param stat pre-calculated statistics that ranks all transcripts
##' after their association with the metabolites (genes in rows)
##' @param members numerical indicating the genes that should be plotted
##' @param toprat the top ratio of genes that should receive the
##' highest color intensity
##' @param rows display layout, rows
##' @param cols display layout,k columns
##' @param cex.text expansion factor for the text
##' @param ... passed on to \code{par}
##' @examples
##' data(tdat)
##' data(mdat)
##' library(org.At.tair.db)
##' pearson <- assocStat(tdat, mdat, absolute=FALSE)
##' members <- na.omit(match(get("00052", org.At.tairPATH2TAIR), fData(tdat)$agi))
##' colorLineplot(tdat, mdat, pearson, members, rows=2, cols=3)
##' @export
##' @return nothing
##' @author Henning Redestig
colorLineplot <- function(tdat, mdat, stat, members, toprat=0.05,
                          rows=NULL, cols=NULL, cex.text=1.1, ...) {
  if(dim(mdat)[1] != dim(stat)[2])
    stop("need one column of statistics for every metabolite")
  x <- tdat$time
  if(is.null(x))
    stop("no time variable in the transcripts data")
  
  fixCol <- function(x, alpha) {
    res <- rep(NA, length(x))
    for(i in 1:length(x)) {
      a <- col2rgb(x[i]) / 255
      res[i] <- rgb(a[1], a[2], a[3], alpha[i])
    }
    names(res) <- names(x)
    res
  }


  rnk <- apply(-stat, 2, rank, ties.method="max")
  par(mfrow=c(ifelse(is.null(rows), length(methods), rows),
        ifelse(is.null(cols), nrow(mdat), cols)), xpd=TRUE, mar=c(5,4,4,5),
      ...)
  gn <- t(scale(t(exprs(tdat[members,]))))
  mn <- t(scale(t(exprs(mdat))))
  for(i in 1:nrow(mdat)) {
    we <- c(1 / rnk[,i] * (nrow(tdat) * toprat))
    we[we > 1] <- 1
    mmat <- rbind(gn,mn[i,])
    rnkm <- c(rnk[members,i], NA)
    allColors <- fixCol(c(rainbow(length(members)),"#000000"), c(we[members],1))
    names(allColors) <- c(rownames(gn), makeMetLabels2(mdat)[i])
    matplot(replicate(nrow(mmat), x),
            y=t(mmat),lwd=2,
            type="l",lty=1,
            col=allColors,
            ylab="Expression/Abundance",
            xlab="Time (h)", xaxt="n", xlim=c(min(x), max(x) * 1.25),
            bty="n")
    delta <- diff(range(gn, na.rm=TRUE)) / nrow(gn)
    o <- order(mmat[,ncol(mmat)])
    axis(1, at=x)
    floor <- function(x, y) { x[ x > y ] <- y; x}
    for(j in 1:length(o)) {
      lines(c(max(x) * 1.007, max(x) * 1.01, max(x) * 1.05, max(x) * 1.06),
            c(mmat[o[j],ncol(mmat)],
              mmat[o[j],ncol(mmat)],
              min(mmat, na.rm=TRUE) + delta * j,
              min(mmat, na.rm=TRUE) + delta * j),
            col=allColors[o[j]])
      text(max(x) * 1.06, min(mmat, na.rm=TRUE) + delta * j,
           paste(paste(names(allColors)[o[j]],
                       ifelse(is.na(rnkm[o[j]]),"", rnkm[o[j]]),
                       sep=ifelse(is.na(rnkm[o[j]]),"", ", "))),
           pos=4, cex=cex.text,
           col=allColors[o[j]])
    }
  }
}

