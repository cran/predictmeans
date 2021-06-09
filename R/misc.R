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
  if (inherits(x, "dist")) 
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
  LetMat <- array(TRUE, dim = c(n, 1), dimnames = list(Lvls, NULL))
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
# https://rpubs.com/bbolker/waldvar

# waldVar2 <- function(object) {
  # ## test for/warn if ML fit?
  # dd <- lme4::devfun2(object,useSc=TRUE,signames=FALSE)
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
######################
# vcov.VarCorr.merMod <- function(object,fit,...) {
  # if (isREML(fit)) {
    # warning("refitting model with ML")
    # fit <- refitML(fit)
  # }
  # if (!require("numDeriv")) stop("numDeriv package required")
  # useSc <- attr(object,"useSc")
  # dd <- lme4:::devfun2(fit,useSc=useSc,signames=FALSE)
  # vdd <- as.data.frame(object,order="lower.tri")
  # pars <- vdd[,"sdcor"]
  # npar0 <- length(pars)
  # if (isGLMM(fit)) {
    # pars <- c(pars,fixef(fit))
  # }
  # hh1 <- hessian(dd,pars)
  # vv2 <- 2*solve(hh1)
  # if (isGLMM(fit)) {
    # vv2 <- vv2[1:npar0,1:npar0,drop=FALSE]
  # }
  # nms <- apply(vdd[,1:3],1,
               # function(x) paste(na.omit(x),collapse="."))
  # dimnames(vv2) <- list(nms,nms)
  # return(vv2)
# }

# http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html

#######################################################
# from package merDeriv vcov.lmerMod.R

vcov_lmerMod <- function (object, ...) {  
    if (!is(object, "lmerMod")) 
        stop("vcov.lmerMod() only works for lmer() models.")
    dotdotdot <- list(...)
    if ("full" %in% names(dotdotdot)) {
        full <- dotdotdot$full
    }
    else {
        full <- FALSE
    }
    if ("information" %in% names(dotdotdot)) {
        information <- dotdotdot$information
    }
    else {
        information <- "expected"
    }
    if (!(full %in% c("TRUE", "FALSE"))) 
        stop("invalid 'full' argument supplied")
    if (!(information %in% c("expected", "observed"))) 
        stop("invalid 'information' argument supplied")
    if ("ranpar" %in% names(dotdotdot)) {
        ranpar <- dotdotdot$ranpar
    }
    else {
        ranpar <- "var"
    }
    parts <- getME(object, "ALL")
    yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
    uluti <- length(parts$theta)
    Zlam <- tcrossprod(parts$Z, parts$Lambdat)
    V <- (tcrossprod(Zlam, Zlam) + Matrix::Diagonal(parts$n, 1)) * (parts$sigma)^2
    M <- solve(chol(V))
    invV <- tcrossprod(M, M)
    LambdaInd <- parts$Lambda
    LambdaInd@x[] <- parts$Lind
    invVX <- crossprod(parts$X, invV)
    Pmid <- solve(crossprod(parts$X, t(invVX)))
    P <- invV - tcrossprod(crossprod(invVX, Pmid), t(invVX))
    fixvar <- solve(tcrossprod(crossprod(parts$X, invV), t(parts$X)))
    if (full == FALSE) {
        fixvar
    }
    else {
        fixhes <- tcrossprod(crossprod(parts$X, invV), t(parts$X))
        uluti <- length(parts$theta)
        devV <- vector("list", (uluti + 1))
        devLambda <- vector("list", uluti)
        score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
        for (i in 1:uluti) {
            devLambda[[i]] <- Matrix::forceSymmetric(LambdaInd == i, 
                uplo = "L")
            devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])), 
                parts$Z)
        }
        devV[[(uluti + 1)]] <- Matrix::Diagonal(nrow(parts$X), 1)
        ranhes <- matrix(NA, nrow = (uluti + 1), ncol = (uluti + 
            1))
        entries <- rbind(matrix(rep(1:(uluti + 1), each = 2), 
            (uluti + 1), 2, byrow = TRUE), t(combn((uluti + 1), 
            2)))
        entries <- entries[order(entries[, 1], entries[, 2]), 
            ]
        if (parts$devcomp$dims[["REML"]] == 0) {
            if (information == "expected") {
                ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 
                  1, function(x) as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV, 
                    devV[[x[1]]]), invV), t(devV[[x[2]]])))))
            }
            if (information == "observed") {
                ranhes[lower.tri(ranhes, diag = TRUE)] <- unlist(apply(entries, 
                  1, function(x) as.vector(-as.numeric((1/2) * 
                    lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV, 
                      devV[[x[1]]]), invV), t(devV[[x[2]]])))) + 
                    tcrossprod((tcrossprod((crossprod(yXbe, tcrossprod(tcrossprod(crossprod(invV, 
                      devV[[x[1]]]), invV), t(devV[[x[2]]])))), 
                      invV)), t(yXbe)))))
            }
        }
        if (parts$devcomp$dims[["REML"]] > 0) {
            if (information == "expected") {
                ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 
                  1, function(x) as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(P, 
                    devV[[x[1]]]), P), t(devV[[x[2]]])))))
            }
            if (information == "observed") {
                ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 
                  1, function(x) -as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(P, 
                    devV[[x[1]]]), P), t(devV[[x[2]]])))) + tcrossprod((tcrossprod((crossprod(yXbe, 
                    tcrossprod(tcrossprod(crossprod(invV, devV[[x[1]]]), 
                      P), t(devV[[x[2]]])))), invV)), t(yXbe)))
            }
        }
        ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
        if (information == "expected") {
            varcov_beta <- matrix(0, length(devV), length(parts$beta))
        }
        if (information == "observed") {
            varcov_beta <- matrix(NA, length(devV), length(parts$beta))
            for (j in 1:(length(devV))) {
                varcov_beta[j, ] <- as.vector(tcrossprod(crossprod(parts$X, 
                  (tcrossprod(crossprod(invV, devV[[j]]), invV))), 
                  t(yXbe)))
            }
        }
        if (ranpar == "var") {
            ranhes <- ranhes
            varcov_beta <- varcov_beta
        }
        else if (ranpar == "sd") {
            sdcormat <- as.data.frame(VarCorr(object, comp = "Std.Dev"), 
                order = "lower.tri")
            sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- sdcormat$sdcor[which(is.na(sdcormat$var2))] * 
                2
            sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- sdcormat$vcov[which(!is.na(sdcormat$var2))]/sdcormat$sdcor[which(!is.na(sdcormat$var2))]
            varcov_beta <- sweep(varcov_beta, MARGIN = 1, sdcormat$sdcor2, 
                `*`)
            weight <- apply(entries, 1, function(x) sdcormat$sdcor2[x[1]] * 
                sdcormat$sdcor2[x[2]])
            ranhes[lower.tri(ranhes, diag = TRUE)] <- weight * 
                ranhes[lower.tri(ranhes, diag = TRUE)]
            ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
        }
        else {
            stop("ranpar needs to be var or sd for lmerMod object.")
        }
        full_varcov <- solve(rbind(cbind(fixhes, t(varcov_beta)), 
            cbind(varcov_beta, ranhes)))
        colnames(full_varcov) <- c(names(parts$fixef), paste("cov", 
            names(parts$theta), sep = "_"), "residual")
        callingFun <- try(deparse(sys.call(-2)), silent = TRUE)
        if (length(callingFun) > 1) 
            callingFun <- paste(callingFun, collapse = "")
        if (!inherits(callingFun, "try-error") & grepl("summary.merMod", 
            callingFun)) {
            return(fixvar)
        }
        else {
            return(full_varcov)
        }
    }
}

lav_matrix_trace <- function (..., check = TRUE) 
{
  if (nargs() == 0L) 
    return(as.numeric(NA))
  dots <- list(...)
  if (is.list(dots[[1]])) {
    mlist <- dots[[1]]
  }
  else {
    mlist <- dots
  }
  nMat <- length(mlist)
  if (nMat == 1L) {
    S <- mlist[[1]]
    if (check) {
      stopifnot(NROW(S) == NCOL(S))
    }
    out <- sum(S[lav_matrix_diag_idx(n = NROW(S))])
  }
  else if (nMat == 2L) {
    out <- sum(mlist[[1]] * t(mlist[[2]]))
  }
  else if (nMat == 3L) {
    A <- mlist[[1]]
    B <- mlist[[2]]
    C <- mlist[[3]]
    B2 <- B %*% C
    out <- sum(A * t(B2))
  }
  else {
    M1 <- mlist[[1]]
    M2 <- mlist[[2]]
    for (m in 3L:nMat) {
      M2 <- M2 %*% mlist[[m]]
    }
    out <- sum(M1 * t(M2))
  }
  out
}

lav_matrix_diag_idx <- function (n = 1L) 
{
  1L + (seq_len(n) - 1L) * (n + 1L)
}



