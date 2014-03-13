covariatemeans <- function (model, modelterm, covariate, level=0.05, Df, trans, responsen, plot=TRUE, plotord, mtitle, jitterv=0, newwd=TRUE) {

  if (missing(modelterm)) modelterm <- covariate

  vars <- unlist(strsplit(modelterm, "\\:"))
  ctr.matrix <- Kmatrix(model, modelterm, covariate)
  KK <- ctr.matrix$K
  pltdf <- ctr.matrix$fctnames
  response <- ctr.matrix$response
  mp <- mymodelparm(model)
  bhat <- mp$coef
  
  # We'll work only with the non-NA elements of bhat
  KK <- KK[, mp$estimable, drop=FALSE]   
  pltdf$yhat <- KK%*%bhat
  pltdf$ses <- sqrt(base::diag(KK %*% tcrossprod(mp$vcov, KK)))

  if (missing(Df)) {
  	if (class(model)[1] == "lme") {
          Df <- terms(model$fixDF)[modelterm]
        }else if (class(model)[1] == "lmerMod") {
      		termlabel <- attr(terms(model),"term.labels")
      		for (i in vars) termlabel <- termlabel[grep(i, termlabel)]
      		termlabel <- paste(termlabel, collapse="-")
      		model.b <- update( model, as.formula(paste(".~. -", termlabel)))
      		Df <- getKR(KRmodcomp(model, model.b), "ddf")
        }else Df <- mp$df
  	if (Df==0) stop("You need provide Df for this model!")
   }
  
  if (class(model)[1]=="glm" && missing(trans)) trans <- model$family$linkinv
  if (class(model)[1] == "glmerMod" && missing(trans)) trans <- slot(model, "resp")$family$linkinv

  Mean <- LL <- UL <- xvar <- factors <- bky <- NULL
  if (missing(trans)) {
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses
  }else{
    pltdf$Mean <- trans(pltdf$yhat)
    pltdf$LL <- trans(pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses)
    pltdf$UL <- trans(pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses)
  }
  pltdf$yhat <- pltdf$ses <- NULL 
  if (modelterm==covariate) pltdf$factors <- factor(1) else pltdf$factors <- factor(do.call("paste", c(pltdf[, vars, drop=FALSE], sep=":")))
  colnames(pltdf)[colnames(pltdf)==covariate] <- "xvar"
  
  # delete empty factor combinations  
  mdf <- model.frame(model)
  if (!(response %in% names(mdf))) mdf[,response] <- eval(parse(text=response), mdf)
  
  ndf <- data.frame(table(mdf[, vars, drop = FALSE]))
  if (any(ndf$Freq==0)) { 
    ndf0 <- ndf[ndf$Freq==0, , drop=FALSE] 
    ndf0$factors <- factor(do.call("paste", c(ndf0[, vars, drop=FALSE], sep=":")))
    pltdf <- pltdf[!pltdf$factors%in%ndf0$factors, ]
  }   
  
  if (plot) {
     if (missing(mtitle)) mtitle <- paste("Fitted and observed relationship with", (1-level)*100, "% CI")
     if (missing(trans)) {
       mdf$bky <- mdf[, response]
     }else{
       if (response %in% names(mdf)) {    ## Transformed y before modelling
         if (class(model)[1]%in%c("glm", "glmerMod")) {
           if (class(mdf[, response])=="factor") {
               mdf$bky <- as.numeric(mdf[, response])-1
             }else if (!is.null(dim(mdf[, response]))) {
               mdf$bky <- mdf[, response][,1]/rowSums(mdf[, response])
               if (missing(responsen)) stop("Please provide suitable name for response variable using option 'responsen'!")
               response <- responsen
             }else mdf$bky <- mdf[, response]
          if (isTRUE(all.equal(trans,function(x) x))) {
            if (class(model)[1]=="glm") mdf$bky <- model$family$linkfun(mdf$bky)
            if (class(model)[1] == "glmerMod") mdf$bky <- slot(model, "resp")$family$linkfun(mdf$bky)
            if (missing(responsen)) stop("Please provide suitable name for response variable using option 'responsen'!")
            response <- responsen
          }
         }else{
           mdf$bky <- trans(mdf[, response])
           if (missing(responsen))  stop("Please provide suitable name for response variable using option 'responsen'!")
           response <- responsen
         }
       }else{       ## Transformed y within modelling
         if (missing(responsen))  stop("Please provide suitable name for response variable using option 'responsen'!")
         response <- responsen
         mdf$bky <- mdf[, response]
       }
     }
    if (modelterm==covariate) mdf$factors <- factor(1) else mdf$factors <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
	  names(mdf)[names(mdf)==covariate] <- "xvar"

    if (newwd) dev.new()
    plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), geom="line", ylab=paste(response, "\n"), data=pltdf, main=paste(mtitle, "\n"), colour=factors) +
        geom_smooth(aes(ymin = LL, ymax = UL), data=pltdf, stat="identity")+
        geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)+
    	  guides(col = guide_legend(modelterm))+
    	  theme_bw()
 	  if (modelterm==covariate)  plt <- plt+ theme(legend.position="none")
    print(plt)

    if (length(vars)==2) {
      if (newwd) dev.new()
      if (missing(plotord)) plotord <- 1:2
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), ylab=paste(response, "\n"), main=paste(mtitle, "\n"), data=pltdf,
                geom="line", colour=factor(eval(parse(text = fact1)))) +
        facet_grid(eval(parse(text = paste("~",fact2, sep=""))))+
        geom_smooth(aes(ymin = LL, ymax = UL), data=pltdf, stat="identity") +
        geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)+
        guides(col = guide_legend(fact1))+
    	  theme_bw()
      print(plt)
    }
    if (length(vars)==3) {
      if (newwd) dev.new()
      if (missing(plotord)) plotord <- 1:3
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      fact3 <- (vars[plotord])[3]
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), ylab=paste(response, "\n"), main=paste(mtitle, "\n"), data=pltdf,
                geom="line", colour=factor(eval(parse(text = fact1)))) +
        facet_grid(eval(parse(text = paste(fact2, "~",fact3, sep=""))))+
        geom_smooth(aes(ymin = LL, ymax = UL), data=pltdf, stat="identity") +
        geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)+
        guides(col = guide_legend(fact1))+
    	  theme_bw()
      print(plt)
    }
  }

  return(invisible(plt))
}
