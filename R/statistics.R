##' Compute other (than HMM based) statistics for associations between
##' transcripts and metabolites.
##' @title Other statistics
##' @param tdat transcripts ExpressionSet
##' @param mdat metabolites ExpressionSet
##' @param method the method to use
##' @param absolute logical indicting if inverted patterns also be
##' considered. For ordinary Pearson, this implies concatenating the
##' absolute correlations to the result.
##' @param ... arguments passed on to the chosen method, i.e.,
##'ode{cor} or \code{lagcor} or \code{trainrank} (for the HMM)
##' @return A matrix with statistics
##' @examples
##' data(tdat)
##' data(mdat)
##' library(org.At.tair.db)
##' stat <- assocStat(tdat, mdat, use="pairwise")
##' members <- rep(FALSE, nrow(tdat))
##' members[match(get("00052", org.At.tairPATH2TAIR), fData(tdat)$agi, nomatch=0)] <- TRUE
##' cstat <- consStat(stat, members)
##' boxplot(cstat$cct~members)
##' @export
##' @author Henning Redestig
assocStat <- function(tdat, mdat, method=c("pearson","spearman","lagpearson", "hmm"),
                      absolute=TRUE, ...) {
  method <- match.arg(method)
  mm <- t(exprs(mdat))
  tt <- t(exprs(tdat))
  stat <- switch(method,
                 pearson={
                   tmp <- t(cor(mm, tt, ...))
                   if(absolute)
                     tmp <- cbind(tmp, abs(tmp))
                   tmp
                 },
                 spearman={
                   tmp <- t(cor(apply(mm,2,rank), apply(tt,2,rank), ...))
                   if(absolute)
                     tmp <- cbind(tmp, abs(tmp))
                   tmp
                 },
                 lagpearson={
                   tt[tt == 0] <- NA
                   tmp <- apply(mm, 2, function(x) {
                     apply(tt, 2,
                           function(y) {
                             lagcor(x, y, ...)
                           })
                   })
                   if(absolute)
                     tmp <- cbind(tmp, abs(tmp))
                   tmp
                 },
                 hmm={
                   tr <- trainrank(t(mm), t(tt))
                   tmp <- as.matrix(tr$positive)
                   if(absolute)
                     tmp <- cbind(tmp, as.matrix(tr$negative))
                   tmp
                 })
  if(!absolute) colnames(stat) <- makeMetLabels2(mdat) else
  colnames(stat) <- c(makeMetLabels2(mdat),
                      paste("n", makeMetLabels2(mdat), sep="_"))
  stat
}

##' (best) Lagged pearson
##' @title Lagged pearson
##' @param x a numeric vector
##' @param y a numeric vector
##' @param lags lags to consider for the x-vector
##' @param ... passed on to cor
##' @return the best correlation over the considered lags
##' @examples
##' lagcor(rnorm(10), rnorm(10), -3:3)
##' @export
##' @author henning
lagcor <- function(x, y, lags=-2:2, ...) {
  nx <- length(x)
   xx <- sapply(lags, function(l) {
     if(l == 0) return(x)
     na <- rep(NA, abs(l))
     if(l > 0)
       return(c(na, x[-((nx-(l-1)):nx)])) else
     return(c(x[-(1:abs(l))], na))
   })
  cors <- cor(y, xx,  use="pair", ...)
  tmp <- cors[which.max(abs(cors))]
  ifelse(length(tmp) > 0, tmp, NA)
}

##' Get cross-validation segments that have (as far as possible) the
##' same ratio of all classes.
##' @title Class balanced CV segments
##' @param x a factor, character or numeric vector that describes
##' class membership of a set of items.
##' @param fold the desired number of segments
##' @param seed randomization seed for reproducibility
##' @return a matrix where each row is CV segment
##' @examples
##' seg <- getcvseg(iris$Species, 10)
##' apply(seg, 1, function(s) table(iris$Species[s]))
##' @export
##' @author Henning Redestig
getcvseg <- function(x, fold=7, seed=123) {
  if(any(table(x) < fold)) {
    fold <- min(table(x))
  }
  if(fold < 2)
    stop("too few observations in the smallest class")
  res <-
    sapply(unique(x), function(z) {
      if(!is.null(seed))
        set.seed(seed)
      tmp <- sample(which(x == z))
      seg <- matrix(c(tmp,
                      rep(NA, ifelse(length(tmp) %% fold ==0, 0,
                                     fold - (length(tmp) %% fold)))),
                    nrow=fold)
    },simplify=FALSE)
  res <- do.call("cbind", res)
  res[!apply(is.na(res), 1, all),,drop=FALSE]
}

##' Calculate a consensus statistic from a number of statistics that
##' each measure the association between a set of metabolites and a
##' set of transcripts. 
##'
##' With an matrix of association statistics between a set of
##' metabolites and transcripts this function uses a set of known
##' assocaition to perform a supervised summarization using the
##' PLS+CCA replacement for OPLS suggested by MacGregor. Briefly, a
##' PLS-DA model is calculated between using the known true
##' associations as guiding examples and PLS scores are
##' extracted. These are then rotated using  canonical covariate
##' analysis to obtain the component that is truly correlated to class
##' separation. The resulting score vector is the consensus statistic
##' (which is high for strong associations but has not other
##' interpretation).
##' @title Consensus statistic
##' @param stat a matrix of association statistics (rows are transcripts)
##' @param members a logical vector indiating if the transcript is
##' associated to the metabolite(s).
##' @param cv Do cross-validation and return the cross-predicted
##' @param fold Number of CV-folds
##' @return A list with the consesus statistic (\code{cct}), the
##' correlation loadings (\code{ccp}) and possibly the cv predicted
##' statistic (\code{ccthat}).
##' @export
##' @references Yu, H and MacGregor, JF. (2004) Post processing
##' methods (PLS-CCA): simple alternatives to pre-processing methods
##' (OSC-PLS). Chemometrics and Intelligent Laboratory Systems 73,
##' 199-205.
##' @examples
##' stat <- matrix(rnorm(100*5), ncol=5)
##' members <- c(rep(1, 21), rep(0, 79))
##' consStat(stat, members)
##' @author Henning Redestig
consStat <- function(stat, members, cv=FALSE, fold=7) {
  mm <- model.matrix(~-1+factor(members))
  colnames(mm) <- levels(factor(members))
  mmm <- scale(as.integer(members))
  pfit <- plsr(mmm~stat, ncomp=ncol(stat)-1, validation="CV",
               segments=apply(getcvseg(members, fold), 1, na.omit))
  nc <- guessComp2(pfit, mmm)$best
  cc <- cancor(pfit$scores[,1:nc], mmm)
  cct <- (pfit$scores[,1:nc] %*% cc$xcoef[,1,drop=FALSE]) * drop(ifelse(cc$ycoef > 0, 1, -1))
  ccp <- drop(cor(cct, stat))
  if(cv) {
    cvseg <- getcvseg(members, fold)
    ccthat <- cct * 0
    for(i in 1:nrow(cvseg)) {
      spfit <- plsr(mmm[-na.omit(cvseg[i,])]~stat[-na.omit(cvseg[i,]),], ncomp=nc)
      scc <- cancor(spfit$scores, mmm[-na.omit(cvseg[i,])])
      plsT <- predict(spfit, stat[na.omit(cvseg[i,]),],type="scores")
      ccthat[na.omit(cvseg[i,]),] <-
        (plsT %*% scc$xcoef[,1,drop=FALSE]) * drop(ifelse(scc$ycoef > 0, 1, -1))
    }
    return(list(cct=cct, ccp=ccp, ccthat=ccthat))
  }
  return(list(cct=cct, ccp=ccp))
}

##' Guess optimal number of components
##'
##' Use some heuristics to guess the best number of components for PLS
##' model. Specifically, take the best component yielding the best Q2
##' and then remove components as long as there is no sudden decrease
##' in Q2.
##' @param fit a \code{plsr} model
##' @param y the y data
##' @param lev the ratio of Q2 that is considered small a decrease.
##' @param minq2 the minimum Q2 to consider, less than this and zero is
##' returned instead.
##' @return a list with R2, Q2 and the recommended number of components
##' @export
##' @examples
##' data(oliveoil)
##' sens.pls <- plsr(sensory~chemical, validation="CV", ncomp=5, data=oliveoil)
##' with(oliveoil, guessComp2(sens.pls, sensory))
##' @author Henning Redestig
guessComp2 <- function (fit, y, lev = 0.05, minq2 = 0) {
    q2 <- drop(R2(fit)$val)[-1]
    sy <- scale(y, center = TRUE, scale = !is.null(fit$scale))
    r2 <- sapply(1:fit$ncomp, function(i) 1 - sum(resid(fit)[, , i]^2)/sum(sy^2))
    best <- which.max(q2)
    if (any(is.na(q2))) 
        best <- 0
    else if (max(q2) > minq2) {
        repeat {
            if (best == 1) 
                break
            if ((q2[best] - q2[best - 1])/q2[best] >= lev) 
                break
            best <- best - 1
        }
    }
    else best <- 0
    list(best = best, q2 = q2, r2 = r2)
}

accuracy <- function(hat, tru) {
  tp <- sum(hat & tru)
  tn <- sum(!hat & !tru)
  fp <- sum(hat & !tru)
  fn <- sum(!hat & tru)
  (tp + tn) / (tp + tn + fp + fn)
}

svdi <- function(x, threshold=0.01, maxSteps=100, ...) {
  missing <- is.na(x)
  temp <- apply(missing, 2, sum)
  missIx <- which(temp != 0)
  x[missing] <- 0
  count <- 0
  error <- Inf
  while ((error > threshold) && (count < maxSteps)) {
    res <- prcomp(t(x), center=FALSE, scale=FALSE, retx=TRUE)
    loadings <- res$rotation
    sDev <- res$sdev
    for (index in missIx) {
      target <- x[!missing[, index], index, drop = FALSE]
      Apart <- loadings[!missing[, index], , drop = FALSE]
      Bpart <- loadings[missing[, index], , drop = FALSE]
      X <- ginv(Apart) %*% target
      estimate <- Bpart %*% X
      x[missing[, index], index] <- estimate
    }
    count <- count + 1
    if (count > 5)
      error <- sqrt(sum((xOld - x)^2)/sum(xOld^2))
    xOld <- x
  }
  svd(x, ...)
}

oda <- function(x, y, ncomp=NULL, nocomp=0, scale="uv", cv=FALSE) {
  if(cv) 
    seg <- getcvseg(y, ifelse(is.numeric(cv), cv, 7)) else
  seg <- rbind(1:nrow(x))

  ym <- model.matrix(~-1+y)
  if(is.factor(y))
    colnames(ym) <- levels(ym) else
  colnames(ym) <- unique(ym)
  if(is.null(ncomp)) ncomp <- ncol(ym) - 1
  if(ncomp > ncol(ym)) ncomp <- ncol(ym) - 1

  xp <- prep(x, scale, simple=FALSE)
  yp <- prep(ym, scale, simple=FALSE)

  ssy <- sum(yp$data^2, na.rm=TRUE)
  ssx <- sum(xp$data^2, na.rm=TRUE)

  for(i in 1:nrow(seg)) {
    keep <- (1:nrow(x))[!(1:nrow(x)) %in% na.omit(seg[i,])]
    xi <- x[keep,,drop=FALSE]
    yi <- ym[keep,,drop=FALSE]

    xr <- xis <- prep(xi, "uv")
    yr <- yis <- prep(ym, "uv")

    ssxr <- sum(xr^2, na.rm=TRUE)
    ssyr <- sum(yr^2, na.rm=TRUE)
    sigma <- crossprod(yr, xr)
    sres <- svdi(sigma, nu=min(dim(sigma)), nv=min(dim(sigma)))
    for(ci in 1:ncomp) {
      xr <- xis
      yr <- yis
      w <- sres$v[,1:ci,drop=FALSE]
      cc <- sres$u[,1:ci,drop=FALSE]

      tt <- xr %*% w
      ts <- array(0, nrow=c(dim(tt), nocomp))
      ts[,,1] <- tt
      e <- xr - tcrossprod(tt, w)
      tto <- ppo <- wwo <- ttoPrev <- NULL
      if(nocomp > 0) {
        for(oci <- 1:nocomp) {
          tm <- scores(pca(crosspred(e, tt), nPcs=1))
          wo <- tm
          wo <- wo / sqrt(drop(crossprod(wo)))

          to <- xr %*% wo
          po <- crossprod(xr, to / drop(crossprod(to)))

          xr <- xr - tcrossprod(to, po)
          tt <- xr %*% w
          ts[,,(oci + 1)] <- ts
          ttoPrev <- cbind(ttoPrev, tt)
          tto <- cbind(too, to)
          ppo <- cbind(ppo, po)
          woo <- cbind(wwo, wo)

          ep <- e
          e <- xi - tcrossprod(tt, w) - tcrossprod(to, po)
        }
      }
      yr <- yi
      u <- yr %*% cc
      us <- array(0, dim=c(dim(u), 1))
      us[,,1] <- u
      f <- yr - tcrossprod(u, cc)

      if(nrow(seg) == 1){
        rownames(tt) <- rownames(x)
        rownames(u) <- rownames(ym)
        rownames(w) <- colnames(x)
        rownames(cc) <- colnames(ym)
        if(!is.null(to)){
          rownames(to) <- rownames(x)
          rownames(po) <- colnames(x)
        }
      }
      mod <- list(t=tt, w=w, u=u, c=cc, to=to, po=po,
                  xm=xp$center, xs=xp$scale,
                  ym=yp$center, ys=yp$scale)
      if(nrow(seg) > 1) {
        for(i in 1:nocomp) {
          smod <- list(t=tt, w=w, u=u, c=cc, to=to, po=po,
                       xm=xp$center, xs=xp$scale,
                       ym=yp$center, ys=yp$scale)      
                  
        
