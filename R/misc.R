adiag <- function (..., pad = as.integer(0), do.dimnames = TRUE) # function from package 'magic'
{
    args <- list(...)
    if (length(args) == 1) {
        return(args[[1]])
    }
    if (length(args) > 2) {
        jj <- do.call("Recall", c(args[-1], list(pad = pad)))
        return(do.call("Recall", c(list(args[[1]]), list(jj), 
            list(pad = pad))))
    }
    a <- args[[1]]
    b <- args[[2]]
    if (is.null(b)) {
        return(a)
    }
    if (is.null(dim(a)) & is.null(dim(b))) {
        dim(a) <- rep(1, 2)
        dim(b) <- rep(1, 2)
    }
    if (is.null(dim(a)) & length(a) == 1) {
        dim(a) <- rep(1, length(dim(b)))
    }
    if (is.null(dim(b)) & length(b) == 1) {
        dim(b) <- rep(1, length(dim(a)))
    }
    if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
        stop("a and b must have identical number of dimensions")
    }
    s <- array(pad, dim.a + dim.b)
    s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
    ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) + 
        dim.a[[i]])
    out <- do.call("[<-", c(list(s), ind, list(b)))
    n.a <- dimnames(a)
    n.b <- dimnames(b)
    if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
        dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
        names(dimnames(out)) <- names(n.a)
    }
    return(out)
}

vec2mat2 <- function (x, sep = "-") 
{
    splits <- strsplit(x, sep)
    n.spl <- sapply(splits, length)
    if (any(n.spl != 2)) 
        stop("Names must contain exactly one '", sep, "' each;  instead got ", 
            paste(x, collapse = ", "))
    x2 <- t(as.matrix(as.data.frame(splits)))
    dimnames(x2) <- list(x, NULL)
    x2
}

multcompLetters <- function (x, compare = "<", threshold = 0.05,   # function from package 'multcompView'
  Letters = c(letters, LETTERS, "."), reversed = FALSE) 
{
  x.is <- deparse(substitute(x))
  if (class(x) == "dist") 
    x <- as.matrix(x)
  if (!is.logical(x)) 
    x <- do.call(compare, list(x, threshold))
  dimx <- dim(x)
{
    if ((length(dimx) == 2) && (dimx[1] == dimx[2])) {
      Lvls <- dimnames(x)[[1]]
      if (length(Lvls) != dimx[1]) 
        stop("Names requred for ", x.is)
      else {
        x2. <- t(outer(Lvls, Lvls, paste, sep = ""))
        x2.n <- outer(Lvls, Lvls, function(x1, x2) nchar(x2))
        x2.2 <- x2.[lower.tri(x2.)]
        x2.2n <- x2.n[lower.tri(x2.n)]
        x2a <- substring(x2.2, 1, x2.2n)
        x2b <- substring(x2.2, x2.2n + 1)
        x2 <- cbind(x2a, x2b)
        x <- x[lower.tri(x)]
      }
    }
    else {
      namx <- names(x)
      if (length(namx) != length(x)) 
        stop("Names required for ", x.is)
      x2 <- vec2mat2(namx)
      Lvls <- unique(as.vector(x2))
    }
  }
  n <- length(Lvls)
  LetMat <- array(TRUE, dim = c(n, 1), dimnames = list(Lvls, 
                                                       NULL))
  k2 <- sum(x)
  if (k2 == 0) {
    Ltrs <- rep(Letters[1], n)
    names(Ltrs) <- Lvls
    dimnames(LetMat)[[2]] <- Letters[1]
    return(Ltrs)
  }
  distinct.pairs <- x2[x, , drop = FALSE]
  absorb <- function(A.) {
    k. <- dim(A.)[2]
    if (k. > 1) {
      for (i. in 1:(k. - 1)) for (j. in (i. + 1):k.) {
        if (all(A.[A.[, j.], i.])) {
          A. <- A.[, -j., drop = FALSE]
          return(absorb(A.))
        }
        else {
          if (all(A.[A.[, i.], j.])) {
            A. <- A.[, -i., drop = FALSE]
            return(absorb(A.))
          }
        }
      }
    }
    A.
  }
  for (i in 1:k2) {
    dpi <- distinct.pairs[i, ]
    ijCols <- (LetMat[dpi[1], ] & LetMat[dpi[2], ])
    if (any(ijCols)) {
      A1 <- LetMat[, ijCols, drop = FALSE]
      A1[dpi[1], ] <- FALSE
      LetMat[dpi[2], ijCols] <- FALSE
      LetMat <- cbind(LetMat, A1)
      LetMat <- absorb(LetMat)
    }
  }
  sortCols <- function(B) {
    firstRow <- apply(B, 2, function(x) which(x)[1])
    B <- B[, order(firstRow)]
    firstRow <- apply(B, 2, function(x) which(x)[1])
    reps <- (diff(firstRow) == 0)
    if (any(reps)) {
      nrep <- table(which(reps))
      irep <- as.numeric(names(nrep))
      k <- dim(B)[1]
      for (i in irep) {
        i. <- i:(i + nrep[as.character(i)])
        j. <- (firstRow[i] + 1):k
        B[j., i.] <- sortCols(B[j., i., drop = FALSE])
      }
    }
    B
  }
  LetMat. <- sortCols(LetMat)
  if (reversed) 
    LetMat. <- LetMat.[, rev(1:ncol(LetMat.))]
  k.ltrs <- dim(LetMat.)[2]
  makeLtrs <- function(kl, ltrs = Letters) {
    kL <- length(ltrs)
    if (kl < kL) 
      return(ltrs[1:kl])
    ltrecurse <- c(paste(ltrs[kL], ltrs[-kL], sep = ""), 
                   ltrs[kL])
    c(ltrs[-kL], makeLtrs(kl - kL + 1, ltrecurse))
  }
  Ltrs <- makeLtrs(k.ltrs, Letters)
  dimnames(LetMat.)[[2]] <- Ltrs
  LetVec <- rep(NA, n)
  names(LetVec) <- Lvls
  for (i in 1:n) LetVec[i] <- paste(Ltrs[LetMat.[i, ]], collapse = "")
  nch.L <- nchar(Ltrs)
  blk.L <- rep(NA, k.ltrs)
  for (i in 1:k.ltrs) blk.L[i] <- paste(rep(" ", nch.L[i]), 
                                        collapse = "")
  monoVec <- rep(NA, n)
  names(monoVec) <- Lvls
  for (j in 1:n) {
    ch2 <- blk.L
    if (any(LetMat.[j, ])) 
      ch2[LetMat.[j, ]] <- Ltrs[LetMat.[j, ]]
    monoVec[j] <- paste(ch2, collapse = "")
  }
  return(monoVec)
}

######################## Functions from lmerTest #################

calcKRDDF <- function(model, term){
  
  rho <- list() ## environment containing info about model
  rho$model <- model
  rho$fixEffs <- fixef(rho$model)
  rho$sigma <- sigma(rho$model)
  rho$thopt <- getME(rho$model, "theta")
  rho$Xlist <- createDesignMat(rho) ## X design matrix for fixed effects
  
  ## define the terms that are to be tested
  rho$test.terms <- attr(terms(rho$model), "term.labels")[unique(attr(rho$Xlist$X.design.red, "assign"))] 
  ## calculate general set of hypothesis matrix 
  L <- calcGeneralSet12(rho$Xlist$X.design.red)
  L <- makeContrastType1(L, rho$Xlist$X.design.red, 
                         rho$Xlist$trms, rho$test.terms)
  Lc <- L[[term]]			 	
  # DDF <- pbkrtest::get_Lb_ddf(rho$model, Lc)
  res.KR <- pbkrtest::KRmodcomp(rho$model, Lc )
  return(res.KR)
}

calcGeneralSet12 <- function(X){
  p <- ncol(X)
  XtX <- crossprod(X)
  U <- doolittle(XtX)$U
  d <- diag(U)
  for(i in 1:nrow(U))
    if(d[i] > 0) U[i, ] <- U[i, ] / d[i]
  L <- U
  L
}


makeContrastType1 <- function(L, X, trms, trms.lab){
  asgn <- attr(X, "assign")
  #trms.lab <- attr(trms, "term.labels")
  p <- ncol(X)
  ind.list <- split(1L:p, asgn)
  df <- unlist(lapply(ind.list, length))
  Lt.master <- L
  Lt.list <- lapply(ind.list, function(i) Lt.master[i, , drop=FALSE])
  if(attr(trms,"intercept") == 0)
    names(Lt.list) <- trms.lab
  else
    names(Lt.list) <- c("(Intercept)", trms.lab)
  Lt.list
}

########################
## construct design matrix for F test 
createDesignMat <- function(rho)
{
  model.term <- terms(rho$model)
  fixed.term <- attr(model.term,"term.labels") 
  X.design <- names.design.mat <-  names.design <- NULL
  X.design.red <- model.matrix(rho$model)
  attr(X.design.red, "dataClasses") <- 
    attr(terms(rho$model, FALSE), "dataClasses")
  dd <- model.frame(rho$model) 
  
  for(i in 1:length(fixed.term))
  {
    formula.term <- as.formula(paste("~", fixed.term[i], "- 1"))
    X.design <- cbind(X.design, model.matrix(formula.term, dd))
    names.design.mat <- c(names.design.mat, 
                          rep(fixed.term[i],
                              ncol(model.matrix(formula.term, dd))))
  }
  
  if(attr(model.term, "intercept") != 0){
    names.design <- c("(Intercept)", colnames(X.design))
    X.design <- cbind(rep(1,dim(X.design)[1]), X.design)
    names.design.mat <- c("(Intercept)", names.design.mat)
  }
  else
    names.design <- colnames(X.design)
  colnames(X.design) <- names.design.mat
  
  # if(length(which(colSums(X.design)==0)) != 0){
    # warning(paste("missing cells for some factors (combinations of factors) \n", 
                  # "care must be taken with type",
                  # as.roman(3) ,
                  # " hypothesis "))
  # }
  
  fullCoefs <- rep(0, ncol(X.design))
  fullCoefs <- setNames(fullCoefs, names.design) 
  if("(Intercept)" %in% names.design)
    names(fullCoefs)[1] <- "(Intercept)"
  fullCoefs[names(rho$fixEffs)] <- rho$fixEffs
  nums.Coefs <- which(names(fullCoefs) %in% names(rho$fixEffs))
  nums.Coefs <- setNames(nums.Coefs, names(fullCoefs[nums.Coefs])) 
  Xlist <- list(X.design.red = X.design.red, 
                trms = model.term,
                X.design = X.design,
                names.design = names.design,
                fullCoefs = fullCoefs,
                nums.Coefs = nums.Coefs)
  return(Xlist)
}
##################

doolittle <- function(x, eps = 1e-6) {
  
  if(!is.matrix(x)) stop("argument 'x' is not a matrix")
  if(ncol(x) != nrow(x))
    stop( "argument x is not a square matrix" )
  if (!is.numeric(x) )
    stop( "argument x is not numeric" )
  n <- nrow(x)
  L <- U <- matrix(0, nrow=n, ncol=n)
  diag(L) <- rep(1, n)
  for(i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for(j in 1:n) {
      U[i,j] <- x[i,j]
      if (im1 > 0) {
        for(k in 1:im1) {
          U[i,j] <- U[i,j] - L[i,k] * U[k,j]
        }
      }
    }
    if ( ip1 <= n ) {
      for ( j in ip1:n ) {
        L[j,i] <- x[j,i]
        if ( im1 > 0 ) {
          for ( k in 1:im1 ) {
            L[j,i] <- L[j,i] - L[j,k] * U[k,i]
          }
        }
        L[j, i] <- if(abs(U[i, i]) < eps) 0 else L[j,i] / U[i,i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list( L=L, U=U )
}


#######################

# waldVar2 <- function(object) {
  # ## test for/warn if ML fit?
  # dd <- lme4:::devfun2(object,useSc=TRUE,signames=FALSE)
  # nvp <- length(attr(dd,"thopt"))+1 ## variance parameters (+1 for sigma)
  # pars <- attr(dd,"optimum")[seq(nvp)] ## var params come first
  # hh <- numDeriv::hessian(dd,pars)
  # ## factor of 2: deviance -> negative log-likelihood
  # vv <- 2*solve(hh)
  # nn <- tn(object)
  # dimnames(vv) <- list(nn,nn)
  # return(vv)
# }

# tn <- function(object) {
  # c(names(getME(object,"theta")),"sigma")
# }

# confintlmer <- function (object, parm, level = 0.95, ...) 
# {
  # cf <- coef(object)
  # pnames <- names(cf)
  # if (missing(parm)) 
    # parm <- pnames
  # else if (is.numeric(parm)) 
    # parm <- pnames[parm]
  # a <- (1 - level)/2
  # a <- c(a, 1 - a)
  # pct <- format.perc(a, 3)
  # fac <- qnorm(a)
  # ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  # ses <- sqrt(diag(object$vcov))[parm]
  # ci[] <- cf[parm] + ses %o% fac
  # ci
# }

# format.perc <- function (probs, digits) {
  # paste(format(100 * probs, trim = TRUE, scientific = FALSE,
               # digits = digits), "%")
# }

  ###################### for print
print.pdmlist = function(x, ...){
  pos = grep('predictmeansPlot|predictmeansBKPlot|predictmeansBarPlot|p_valueMatrix', names(x))
  x = x[names(x)[-pos]]
  NextMethod()
}




