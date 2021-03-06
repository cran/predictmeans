predictmeans <- function (model, modelterm, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL, 
                          level=0.05, covariate=NULL, letterdecr=TRUE, trans = NULL, transOff=0, responsen=NULL, count=FALSE, 
                          plotord=NULL, plottitle=NULL, plotxlab=NULL, plotylab=NULL, mplot=TRUE, barplot=FALSE, pplot=TRUE, 
                          bkplot=TRUE, plot=TRUE, jitterv=0, basesz=12, prtnum=TRUE, newwd=TRUE, 
                          permlist=NULL, ndecimal=4) 
{  
  options(scipen=6)
  if(any(missing(model), missing(modelterm))) stop("The arguments 'model', and 'modelterm' must be provided!")
  if (!(modelterm %in% attr(terms(model), "term.labels"))) stop(paste("The", modelterm, "must be exactly a term in the model (especially check the order of interaction)."))
  slevel <- level
  predictmeansPlot <- predictmeansBarPlot <- NULL
  if (inherits(model, "aovlist")) stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glm")) {
    trans <- model$family$linkinv  # identical(trans, make.link("log")$linkinv)
    if (model$family$family %in% c("poisson", "quasipoisson")) count=TRUE
  }
  if (inherits(model, "glmerMod")) {
    trans <- slot(model, "resp")$family$linkinv
    #    if (slot(model, "resp")$family$family %in% c("poisson", "quasipoisson", "binomial", "quasibinomial")) count=TRUE
  }
  vars <- unlist(strsplit(modelterm, "\\:"))
  mdf <- model.frame(model)
  if(any(!is.element(vars, names(mdf)[sapply(mdf,is.factor)])))
    stop(sQuote(vars), "all must be factor(s)!")
  
  if (length(vars)==1) atvar <- NULL
  if (!is.null(permlist) && !permlist%in%c("NULL", "")) {pairwise <- TRUE; if (adj=="tukey") stop("The p-value can't be adjusted by Tukey methd!")}                        # option checking
  if (!is.null(atvar) && !atvar%in%c("NULL", "")) pairwise <- TRUE
  if (adj != "none") pairwise <- TRUE
  if (!is.null(plotord) && !plotord%in%c("NULL", "")) plot <- mplot <- TRUE
   
  ctr.matrix <- Kmatrix(model, modelterm, covariate, prtnum=prtnum) 
  KK <- ctr.matrix$K
  label <- ctr.matrix$fctnames
  rownames(label) <- rownames(KK)
  response <- ctr.matrix$response
  mp <- mymodelparm(model) 
  n.table <- table(mdf[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model
  
  K <- KK[, mp$estimable, drop = FALSE]          # To match coef names
  if (any(ndf$Freq==0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rownames(K)%in%rnTrt,]
    label <- label[rownames(label)%in%rnTrt,]
  }
  pm <- K %*% mp$coef
  vcovm <- mp$vcov
  # ses <- sqrt(diag(K %*% tcrossprod(vcovm, K)))
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  mt <- data.frame(pm, ses, label)
  
  # print(mt)
  # names(mt)[-1:-2] <- vars
  # print(mt)
  bkmt <- mt  # for back transformed
  mean.table <- round(xtabs(pm ~ ., mt[, c("pm", vars)], drop.unused.levels = TRUE), ndecimal)
  se.table <- round(xtabs(ses ~ ., mt[, c("ses", vars)], drop.unused.levels = TRUE), ndecimal+1)
  mean.table[!(n.table)] <- NA
  se.table[!(n.table)] <- NA
  if (length(vars) > 1) {
    varsnlevel <- numeric(0)
    for (i in vars) varsnlevel[i] <- nlevels(mdf[, i])
    tbvars <- names(sort(varsnlevel, decreasing = TRUE))
    mean.table <- ftable(mean.table, row.vars =tbvars[1], col.var=tbvars[-1])
    se.table <- ftable(se.table, row.vars =tbvars[1], col.var=tbvars[-1])
  }
  if (length(na.omit(unique(se.table)))==1) {
    se.table <- min(se.table, na.rm=TRUE)
    names(se.table) <- "All means have the same Stder"
  }
  
  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)
  varn1 <- varn2 <- rep(0, nK * (nK - 1)/2)
  
  if (nK == 1) {
    SED.out <- NA
    LSD <- NA
  }else {
    kindx <- 1:nK
    CM <-  matrix(0, nrow=nK * (nK - 1)/2, ncol=nK)
    t <- 1
    for (i in 2:nK) {                      # To construct pairwise comparison K matrix by col order
      for (j in 1:(i-1)) {
        CM[t, ] <- (kindx == j) - (kindx == i)
        varn1[t] <- rnK[i]
        varn2[t] <- rnK[j]
        t <- t+1
      }
    }
    
    nKK <- nrow(KK)          # To setup various matrix and row, col names
    rnKK <- rownames(KK)
    KKvarn1 <- KKvarn2 <- rep(0, nKK * (nKK - 1)/2)
    tt <- 1
    for (i in 2:nKK) {                      
      for (j in 1:(i-1)) {
        KKvarn1[tt] <- rnKK[i]
        KKvarn2[tt] <- rnKK[j]
        tt <- tt+1
      }
    }
    KKvarndiff <- data.frame(matrix(unlist(strsplit(KKvarn1, "\\:")), byrow=T, nrow=length(KKvarn1)),
                          matrix(unlist(strsplit(KKvarn2, "\\:")), byrow=T, nrow=length(KKvarn2)))
    
    rK <- CM%*%K                    # calculate stats
    cm <- rK%*%mp$coef
    #vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    #dses <- sqrt(diag(vcov.contr))
    dses <- as.numeric(apply(rK, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(vcovm, y))}))
    if (adj == "bonferroni") level <- level/length(dses)
    
    SED.out <- c(Max.SED = max(dses), Min.SED = min(dses), Aveg.SED = mean(dses))
    dses.df <- data.frame(matrix(unlist(strsplit(varn1, "\\:")), byrow=T, nrow=length(varn1)),
                          matrix(unlist(strsplit(varn2, "\\:")), byrow=T, nrow=length(varn2)), dses)
    
    if (all(length(vars) > 1, SED.out[1]!=SED.out[2])) {
      dses.m <- matrix(0, nrow=3, ncol=length(vars))
      colnames(dses.m) <- vars
      rownames(dses.m) <- c("Aveg.SED", "Min.SED", "Max.SED")
      for (i in 1:length(vars)) {
        varsindx <- as.character(dses.df[,i])==as.character(dses.df[, length(vars)+i])
        if(any(varsindx))  # in case of A/B
          dses.m[,i] <- summary.default(dses.df[varsindx,"dses"])[c(4,1,6)]  else {dses.m[,i] <- NA; atvar <- NULL}
      }
      attr(SED.out, "For the Same Level of Factor") <- dses.m
    }
    
    if (is.null(permlist) || permlist%in%c("NULL", "")) {
      if (length(Df) == 0) {
        if (inherits(model, "lme")) {
          Df <- terms(model$fixDF)[modelterm]
        }else if (inherits(model, "merMod")) {
		  Df <- pbkrtest::getKR(calcKRDDF(model, modelterm), "ddf")
        }else Df <- mp$df
        
        if (Df==0) stop("You need provide Df for this model!")
      }
      LSD <- round(qt(1 - level/2, df = Df) * SED.out, ndecimal+1)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "For the Same Level of Factor") <- NULL
      #    if (length(vars) > 1) {
      #      rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
      #      attr(LSD, "For the Same Level of Factor") <- qt(1 - level/2, df = Df) * dses.m
      #    } # end of if LSD
      attr(LSD, "Significant level") <- slevel
      attr(LSD, "Degree of freedom") <- round(Df, 2)
    }else{
      LSD <- round(2 * SED.out[1:3], ndecimal+1)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "Note") <- "This is a approximate LSD which is 2*SED."          
    }
    
	
    if (pairwise) {
	  p_valueMatrix <- NULL
      tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
      rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
      t.v <- cm/dses
      if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj=="tukey")) p.tukey <- ptukey(sqrt(2)*abs(t.v), nK, Df, lower.tail=FALSE)
      tvm[upper.tri(tvm)] <- t.v
      if (is.null(permlist) || permlist%in%c("NULL", "")) {
        t.p.values <- 2 * pt(-abs(t.v), Df)
      }else{      
        nsim <- length(permlist[[1]])
        tValue <- function(x, rK){
          cm <- rK%*%x$coef
          vcovm <- x$vcov
          #vcov.contr <- rK %*% tcrossprod(vcovm, rK)
          #ses <- sqrt(diag(vcov.contr))
          ses <- as.numeric(apply(rK, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(vcovm, y))}))
          t.v <- cm/ses
          return(t.v)	
        }
        if (nK==2)	t.p.values <-  (sum(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)))+1)/(nsim+1) else	
          t.p.values <-  (rowSums(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)))+1)/(nsim+1)
      } # end of if (is.null(permlist))
      
      if (is.null(atvar) || atvar%in%c("NULL", "")) {
        if (adj=="tukey") t.p.valuem[upper.tri(t.p.valuem)] <- p.tukey else
          t.p.valuem[upper.tri(t.p.valuem)] <- p.adjust(t.p.values, adj)
        t.p.valuep <- t.p.valuem    # for plot
        t.p.valuem <- t(t.p.valuem) + tvm
        names(t.p.valuem) <- NULL
        if (is.null(permlist) || permlist%in%c("NULL", "")) {
          attr(t.p.valuem, "Degree of freedom") <- Df
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, p-value (adjusted by",
                                            sQuote(adj), "method) below the diagonal")
          LSDm.up <- qt(1-level/2, df = Df)*dses
          LSDm[upper.tri(LSDm)] <- LSDm.up
          Diffm[upper.tri(Diffm)] <- cm
          LSDm <- t(LSDm)+Diffm
          names(LSDm) <- NULL
          attr(LSDm,"Significant level") <- slevel
          attr(LSDm,"Degree of freedom") <- Df
          attr(LSDm,"Note") <- paste("LSDs matrix has mean differences (row-col) above the diagonal, LSDs (adjusted by",
                                     sQuote(adj), "method) below the diagonal")
        }else{
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, and", sQuote(nsim), "times permutation p-value (adjusted by",
                                            sQuote(adj), "method) below the diagonal")       
        } # end of if (is.null(permlist)) 
		
		groupRn <- rnK[order(mt$pm, decreasing = letterdecr)]
		t.p.valuemGrp <- t.p.valuem
		t.p.valuemGrp[upper.tri(t.p.valuemGrp)] <- t(t.p.valuemGrp)[upper.tri(t.p.valuemGrp)]
		t.p.valuemGrp <- t.p.valuemGrp[groupRn, groupRn]
        # multGrp <- data.frame(multcompLetters(t.p.valuem, Letters=LETTERS, threshold=slevel))
		multGrp <- data.frame(multcompLetters(t.p.valuemGrp, Letters=LETTERS, threshold=slevel)[rnK])
		names(multGrp) <- "Group"
		multGrp <- data.frame(Treatment=groupRn, Mean=mt$pm[order(mt$pm, decreasing = letterdecr)], Group=multcompLetters(t.p.valuemGrp, Letters=LETTERS, threshold=slevel)) 
		rownames(multGrp) <- NULL
         
        attr(t.p.valuem, paste("Letter-based representation of pairwise comparisons at significant level", sQuote(slevel))) <- multGrp
        if (all(nrow(t.p.valuep) > 3, pplot, plot)) {
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Level Plot of p-value (adjusted by", sQuote(adj), "method)\n for Pairwise Comparison")		
          PMplot(t(t.p.valuep), level=slevel, legendx=0.69, mtitle=mtitle, newwd=newwd)
		  p_valueMatrix <- t(t.p.valuep)
        } # end of if (all(nrow(t.p.valuep) > 3, pplot, plot)) 
      }else{
        # atvar <- vars[which(vars %in% atvar)] # ensure a proper order for atvar
        dses.df$tvalue <- t.v
        dses.df$pvalue <- t.p.values        
        dses.df <- merge(dses.df, KKvarndiff, all=TRUE, sort=FALSE)
        atvar.df <- dses.df  
        for (i in which(vars%in%atvar)) { # To find rows relating atvar
          atvar.df <- atvar.df[which(as.character(atvar.df[, i]) == as.character(atvar.df[, length(vars) + i])), ]
        }
        if (adj=="tukey") atvar.df$adj.pvalue <- ptukey(sqrt(2)*abs(atvar.df$tvalue), nK, Df, lower.tail=FALSE) else
          atvar.df$adj.pvalue <- p.adjust(atvar.df$pvalue, adj)
        names(atvar.df)[1:length(vars)] <- vars
        for (i in vars) {                     # To ensure factor vars  have the same level as original
          atvar.df[,i] <- factor(atvar.df[,i], levels=levels(mdf[, i][[1]]))
        }
        atvar.df <- atvar.df[do.call(order, atvar.df[, atvar, drop=FALSE]),]  # sort data by atvar order
        
        rnK.df <- as.data.frame(matrix(unlist(strsplit(rnKK, "\\:")), byrow=T, nrow=nKK)) # To find the suitable names for pm
        colnames(rnK.df) <- vars
        for (i in vars) {       # To ensure factor vars  have the same level as original
          rnK.df[,i] <- factor(rnK.df[,i], levels=levels(mdf[, i][[1]]))
        }
        rnK.df  <- rnK.df [do.call(order, rnK.df[, atvar, drop=FALSE]),]
        resvar <- vars[!vars%in%atvar]   # The rest of vars rather than at var
        rnK.df <- rnK.df[, c(atvar, resvar)]   # To ensure the right matrix name later
        atvar.levels <- unique(do.call("paste", c(rnK.df[, atvar, drop=FALSE], sep=" : ")))
        resvar.levels <- unique(do.call("paste", c(rnK.df[, resvar, drop=FALSE], sep=" : "))) 
        rcnplotm <- do.call("paste", c(rnK.df[, , drop=FALSE], sep=" : ")) # row col names of image plot 
        
		# find decreasing order for each level of atvar mean
        mt.atvar <- mt[do.call(order, mt[, atvar, drop=FALSE]),]
		mt.atvar <- mt.atvar[, c(atvar, resvar, "pm")]
		mt.atvar.factor <- factor(do.call("paste", c(mt.atvar[, atvar, drop=FALSE], sep=" : ")))
		mtlist <- lapply(split(mt.atvar$pm, mt.atvar.factor), function(x) order(x, decreasing = letterdecr))
		mtlist <- mtlist[atvar.levels]
		
        listlength <- length(atvar.levels)
        pmlist <- pmlistTab <- vector("list", listlength)
        indexlength <- nrow(atvar.df)/listlength
        nrow.pm <- length(resvar.levels)
        
        for (i in 1:listlength) {           # extract pvalues for each factor level
          atvar.pm <- matrix(0, nrow=nrow.pm, ncol=nrow.pm)
          atvar.pm[upper.tri(atvar.pm)] <- atvar.df$adj.pvalue[(indexlength*i-(indexlength-1)):(indexlength*i)]
          rownames(atvar.pm) <- colnames(atvar.pm) <- resvar.levels
          pmlist[[i]] <- t(atvar.pm)
		  outtab <- round(as.table(t(atvar.pm)), 4)
		  outtab[outtab < 0.0001] <- "0.0001"
          outtab[col(outtab)==row(outtab)] <- 1.0000
          outtab[upper.tri(outtab)] <-""
		  
		  Grpatvar.pm <- t(atvar.pm)+ atvar.pm
		  Grpatvar.pm <- Grpatvar.pm[mtlist[[i]], mtlist[[i]]]
		  if (all(!is.na(outtab))) outtab <- as.table(cbind(outtab, Group=multcompLetters(Grpatvar.pm, Letters=LETTERS, threshold=slevel)[resvar.levels]))
		  pmlistTab[[i]] <- outtab  		  
        }         
		
        if (all(nrow.pm > 2, pplot, plot)) { 
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Adjusted p-value (by", sQuote(adj), 
                                                                                 "method)\n for Pairwise Comparison at Each Level of",paste(sQuote(atvar), collapse =" and "))
          PMplot(pmlist, level=slevel, xylabel=rcnplotm, legendx=0.69, mtitle=mtitle, newwd=newwd)
		  p_valueMatrix <- pmlist
        }    
      } # if (is.null(atvar))
    }# end of if(pairwise)    

	meanTable <- mt
	meanTable$Df <- Df
	    if (is.null(permlist) || permlist%in%c("NULL", "")) {    
      meanTable$LL <- meanTable$pm - qt(1 - slevel/2, df = Df) * meanTable$ses
      meanTable$UL <- meanTable$pm + qt(1 - slevel/2, df = Df) * meanTable$ses
    }else{
      meanTable$LL <- meanTable$pm - 2 * meanTable$ses
      meanTable$UL <- meanTable$pm + 2 * meanTable$ses
	  slevel <- 0.05
	  meanTable$Df <- NA
    }
	
	rownames(meanTable) <- NULL
	meanTable <- meanTable[c(vars, "pm", "ses", "Df", "LL", "UL")]	
	names(meanTable) <- c(vars, "Predicted means", "Standard error", "Df", paste("LL of ", (1 - slevel) * 100, "% CI", sep = ""),
                                    paste("UL of ", (1 - slevel) * 100, "% CI", sep = ""))
	
    if (plot) {
      if (length(vars) > 3)
        cat("\n", "There is no plot for more than three-way interaction! \n\n")
      plotmt <- na.omit(mt)	  
      yMin <- min(plotmt[, "pm"])
      yMax <- max(plotmt[, "pm"])
      offSet <- 0.25 * (yMax - yMin)
      up <- yMin + LSD[3]
      lsdBar <- cbind(plotmt[, vars, drop=FALSE], up=up, yMin=yMin)
      limits <- aes(ymax = (pm + ses)*(pm > 0) + pmin(pm + ses, 0)*(pm <= 0), ymin=(pm - ses)*(pm < 0) + pmax(pm - ses, 0)*(pm >= 0))
      
      if (length(vars) == 1) {
          mxlab <- plotxlab
          if (is.null(plotxlab) || plotxlab%in%c("NULL", "")) mxlab <- paste("\n", vars, sep="")
          mylab <- plotylab
          if (is.null(plotylab) || plotylab%in%c("NULL", "")) mylab <- paste(response, "\n", sep="")		  
        if (mplot) {
          if (newwd) dev.new()
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for \"", vars, "\" with Aveg.LSD (", slevel * 100, "%) Bar", sep="") 
          p <- qplot(eval(parse(text = vars)), pm, data=plotmt, xlab = mxlab, group=1,
                     ylab = mylab, main = paste(mtitle, "\n", sep=""), 
                     ylim = c(yMin - offSet, max(yMax + offSet, yMin + LSD[3] + offSet)),
                     xlim = c("Aveg.LSD", levels(plotmt[, vars])))
          p <- p +geom_point(colour="red")+
		  geom_line(size=0.5)+
		  geom_errorbar(aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue") +  #data=lsdBar, 
            theme_bw(basesz)
		  predictmeansPlot <- p
          print(p)		  
        } # end of if mplot	   
        if (barplot) {
          if (newwd) dev.new()
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for \"", modelterm, "\" with Stder Bars", sep = "")
          p <- qplot(eval(parse(text = vars)), pm, data=plotmt, ylim = c(min(0, (min(pm) - max(ses)) * 1.2),
                                                                         max(0, (max(pm) + max(ses)) * 1.2)), xlab=mxlab, ylab=mylab, 
                     main=paste(mtitle, "\n", sep=""))
          dodge <- position_dodge(width=0.9)
          p <- p + geom_bar(position=dodge, stat="identity", fill="papayawhip", colour="darkgreen") +
            geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
		  predictmeansBarPlot <- p
          print(p)
        }
      }
      if (length(vars) == 2) {
        if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:2
        fact1 <- (vars[plotord])[1]
        fact2 <- (vars[plotord])[2]
          mxlab <- plotxlab
          if (is.null(plotxlab) || plotxlab%in%c("NULL", "")) mxlab <- paste("\n", fact1, sep="")
          mylab <- plotylab
          if (is.null(plotylab) || plotylab%in%c("NULL", "")) mylab <- paste(response, "\n", sep="")
		
        if (mplot) {     
          if (newwd) dev.new()
          #		if (is.null(plotord)) plotord <- 1:2
          #  fact1 <- (vars[plotord])[1]
          #  fact2 <- (vars[plotord])[2]
          #  nlvel1 <- nlevels(plotmt[, fact1])
          #  nlvel2 <- nlevels(plotmt[, fact2])
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Predicted means for \"", fact1, "\" by \"", fact2, "\" with Aveg.LSD (", 
                                                                                 slevel * 100, "%) Bar", sep = "")
          plotmt[, fact1] <- factor(plotmt[, fact1], levels = c("Aveg.LSD", levels(plotmt[, fact1])))
          p <- qplot(eval(parse(text = fact1)), pm, data=plotmt, group=eval(parse(text = fact2)), 
                     col=eval(parse(text = fact2)), main=paste(mtitle, "\n", sep=""), xlab = mxlab,
                     ylab = mylab, xlim = levels(plotmt[, fact1]),
                     ylim = c(yMin - offSet, max(yMax + offSet, yMin + LSD[3] + offSet)))
          p <- p+ geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), size=0.96)+
            geom_errorbar(aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue")+
            guides(linetype = guide_legend(title = fact2))+
            guides(col = guide_legend(title = fact2))+
            theme_bw(basesz)  
          predictmeansPlot <- p
          print(p)
          #                  theme(legend.position = c(0.12-0.01*(max(nlvel1-5, 0)), 0.88-0.02*(max(nlvel2-3,0))), 
          #                    legend.background = element_rect(colour = "black")) )
        } # end if mplot
        if (barplot) {
          if (newwd) dev.new()
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Predicted means for \"", modelterm, "\" with Stder Bars", sep = "")
          dodge <- position_dodge(width=0.9)
          # p <- qplot(eval(parse(text = fact1)), pm, data=plotmt,
                     # fill=eval(parse(text = fact2)), ylim = c(min(0, (min(pm) - max(ses)) * 1.1),
                                                                              # max(0, (max(pm) + max(ses)) * 1.1)), xlab = mxlab,
                     # ylab = mylab, main = mtitle)+
					 # geom_point(position=dodge) +
					 # geom_bar(stat = "identity", position=dodge)
		  p <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), fill=  eval(parse(text = fact2))))+
            geom_point(position=dodge) +
            geom_bar(stat = "identity", position=dodge)+
			labs( x = mxlab, y = mylab, title =mtitle)+
			ylim(c(min(0, (min(pm) - max(ses)) * 1.1), max(0, (max(pm) + max(ses)) * 1.1)))	 
					 
					 
          p <- p + geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            scale_fill_brewer()+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+            
            theme(legend.position = "top")+guides(fill = guide_legend(title = fact2))
		   predictmeansBarPlot <- p
          print(p)
        }
      }
      
      if (all(length(vars) == 3, mplot)) {
        if (newwd) dev.new()
        mtitle <- plottitle
        if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:3
        fact1 <- (vars[plotord])[1]
        fact2 <- (vars[plotord])[2]
        fact3 <- (vars[plotord])[3]
        plotmt[, fact1] <- factor(plotmt[, fact1], levels = c("Aveg.LSD", levels(plotmt[, fact1])))
        if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for '", fact1, "' by '", fact2, "' for each '",
                                                                              fact3, "'\n with Aveg.LSD (", slevel * 100, "%) Bar\n", sep = "")
          mxlab <- plotxlab
          if (is.null(plotxlab) || plotxlab%in%c("NULL", "")) mxlab <- paste("\n", fact1, sep="")
          mylab <- plotylab
          if (is.null(plotylab) || plotylab%in%c("NULL", "")) mylab <- paste(response, "\n", sep="")
		
																			  
        p <- qplot(eval(parse(text = fact1)), pm,  data=plotmt, main=mtitle,
                   xlab = mxlab, ylab = mylab,
                   xlim = levels(plotmt[, fact1]), ylim=c(yMin-0.5*offSet, max(yMax+0.5*offSet, yMin+LSD[3]+0.5*offSet)),
                   group=factor(eval(parse(text = fact2))), col=factor(eval(parse(text = fact2)))) +
          geom_errorbar(aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue")+
          geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), size=0.8)+
          facet_grid(eval(parse(text = paste("~",fact3, sep=""))))+
          guides(group = guide_legend(fact2))+
          guides(linetype = guide_legend(fact2))+
          guides(col = guide_legend(fact2))+
          theme_bw(basesz)
		 predictmeansPlot <- p
        print(p)
      }
      
    }
  }
  
  if (!is.null(trans)) {
    
    Mean <- Trt <- predictmeansBKPlot <- NULL
	# if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) {
	  # bkmtlog <- t(mapply(bt.log, meanlog=bkmt$pm, sdlog=bkmt$ses, n=rep(round(Df/nrow(bkmt), 0), length(bkmt$pm)), alpha=slevel, SIMPLIFY = T))[, c(1,6,8,9)]
      # bkmtlog[, c(1,3,4)] <- bkmtlog[, c(1,3,4)]-transOff
	  # bkmt$Mean <- bkmtlog[, "btmean"]
	  # #bkmt$Stder <- bkmtlog[, "sd.mean"]
	  # bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = Df) * bkmt$ses)-transOff
	  # bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = Df) * bkmt$ses)-transOff

    # bkmt$pm <- bkmt$ses <- NULL
    # nc <- ncol(bkmt)
    # names(bkmt)[c(nc - 1, nc)] <- c(paste("LL of ", (1 - slevel) * 100, "% CI", sep = ""),
                                    # paste("UL of ", (1 - slevel) * 100, "% CI", sep = ""))
    # bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 4)
    # if (count) bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)
	
	# }else{
	
	
    # bkmt$Mean <- trans(bkmt$pm)
    # if (is.null(permlist) || permlist%in%c("NULL", "")) {    
      # bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = Df) * bkmt$ses)
      # bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = Df) * bkmt$ses)
    # }else{
      # bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)
      # bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)
    # }
	
    # bkmt$pm <- bkmt$ses <- NULL
    # nc <- ncol(bkmt)
    # names(bkmt)[c(nc - 1, nc)] <- c(paste("LL of ", (1 - slevel) * 100, "% CI", sep = ""),
                                    # paste("UL of ", (1 - slevel) * 100, "% CI", sep = ""))
    # bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 4)
    # if (count) bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)
    # }
	
	bkmt$Mean <- trans(bkmt$pm)-transOff
	# if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) bkmt$Mean <- exp(bkmt$pm+bkmt$ses/2)-transOff
	if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) bkmt$Mean <- exp(bkmt$pm)-transOff
    if (is.null(permlist) || permlist%in%c("NULL", "")) {    
      bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = Df) * bkmt$ses)-transOff
      bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = Df) * bkmt$ses)-transOff
    }else{
      bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)-transOff
      bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)-transOff
    }
	
    bkmt$pm <- bkmt$ses <- NULL
    nc <- ncol(bkmt)
    names(bkmt)[c(nc - 1, nc)] <- c(paste("LL of ", (1 - slevel) * 100, "% CI", sep = ""),
                                    paste("UL of ", (1 - slevel) * 100, "% CI", sep = ""))
    bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], ndecimal)
    if (count) {
	  bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)
	  bkmt[, (nc - 2):nc][bkmt[, (nc - 2):nc] < 0] <- 0
	  }	
	
    if (plot && bkplot) {
      if (response %in% names(mdf)) {    ## Transformed y before modelling
        if (inherits(mdf[, response], "factor")){
          bky <- as.numeric(mdf[, response])-1
        }else{
          if (inherits(model, "glm") || inherits(model, "glmerMod")) {
            bky <- mdf[, response]
            if (!is.null(dim(mdf[, response]))) bky <- mdf[, response][,1]/rowSums(mdf[, response])
          }else{
            bky <- trans(mdf[, response])
          }# end of if glm or glmer
        }# end of if factor
      }else{       ## Transformed y within modelling
        nresponse <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
        if (!nresponse %in% names(mdf)) {
          # ques <- paste("R thinks the back transformed response is", sQuote(nresponse), "but the right one is: ")
          # nresponse <- readline(ques)
          if (is.null(responsen) || responsen%in%c("NULL", ""))  stop("Please provide suitable name for response variable using option 'responsen'!")
          nresponse <- responsen
        }
        bky <- mdf[, nresponse]
      }
      if (is.list(bky)) bky <- bky[[1]]
      Trtn <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
      newdata2 <- data.frame(bky=bky, Trtn=Trtn)
      bkmt$Trt <- do.call("paste", c(bkmt[, vars, drop=FALSE], sep=":"))
      xMax <- max(max(bkmt[, nc], na.rm=TRUE), bky, na.rm=TRUE)
      xMin <- min(min(bkmt[, nc - 1], na.rm=TRUE), bky, na.rm=TRUE)
      xoffSet <- 0.15 * (xMax - xMin)
      mtitle <- plottitle
      if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Back Transformed Means with ", (1 - slevel) * 100, "% CIs\n for '",
                                                                            modelterm, "'", "\n", sep = "")
      if (newwd) dev.new()
      limits <- aes(xmax = bkmt$`UL of 95% CI`, xmin=bkmt$`LL of 95% CI`)
      xlimv <- c(xMin - xoffSet, xMax + xoffSet)
      p <- qplot(Mean, Trt, main = mtitle, ylab = "", xlab = "", xlim = xlimv, data=bkmt)+
        geom_point(colour="red") + geom_errorbarh(limits, height=0.2, size=0.8, colour="red") +
        scale_y_discrete(limits = rev(unique(bkmt$Trt)))+
        geom_point(aes(x=bky, y=Trtn), shape=1, position = position_jitter(width = jitterv, height = jitterv), colour="blue", alpha=0.6, data=newdata2)+
        theme_bw(basesz)
	   predictmeansBKPlot <- p
      print(p)
    }
    rownames(bkmt) <- NULL
    bkmt$Trt <- NULL
	
	meanTable <- cbind(meanTable, bkmt[, (ncol(bkmt)-2):ncol(bkmt)])
	colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("BK mean", paste("LL of ", (1 - slevel) * 100, "% BK CI", sep = ""),
                                    paste("UL of ", (1 - slevel) * 100, "% BK CI", sep = ""))
	
	predictmeansPlot <- list(predictmeansPlot, predictmeansBKPlot)
    # if (pairwise) {
      # if (is.null(atvar) || atvar%in%c("NULL", "")) {
        # if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj %in% c("none", "bonferroni"))) {
          # outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      # "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm, ndecimal+1),
                      # "Pairwise p-values" = round(t.p.valuem,4), "Back Transformed Means" = bkmt, 
					  # predictmeansPlot=predictmeansPlot, predictmeansBKPlot=predictmeansBKPlot, 
					  # predictmeansBarPlot=predictmeansBarPlot, p_valueMatrix=p_valueMatrix, mean_table=meanTable)
					  		  # class(outputlist) = "pdmlist"
          # return(outputlist)	
        # }else{
          # outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      # "Standard Error of Differences" = SED.out, LSD = LSD,
                      # "Pairwise p-values" = round(t.p.valuem,4), "Back Transformed Means" = bkmt, 
					  # predictmeansPlot=predictmeansPlot, predictmeansBKPlot=predictmeansBKPlot, 
					  # predictmeansBarPlot=predictmeansBarPlot, p_valueMatrix=p_valueMatrix, mean_table=meanTable)
					  		  # class(outputlist) = "pdmlist"
          # return(outputlist)	
        # }
      # }else{
        # outputlist <- vector("list", listlength+10)
        # outputlist[[1]] <- mean.table
        # outputlist[[2]] <- se.table
        # outputlist[[3]] <- SED.out
        # outputlist[[4]] <- LSD		
		# outputlist[[5]] <- paste("For variable", paste(sQuote(resvar), collapse =" and "), "at each level of", paste(sQuote(atvar), collapse =" and "))                   
		# for (i in 6: (listlength+5))  outputlist[[i]] <- pmlistTab[[i-5]]
        # outputlist[[listlength+6]] <- bkmt
		# outputlist[[listlength+7]] <- predictmeansPlot
		# outputlist[[listlength+8]] <- predictmeansBKPlot
		# outputlist[[listlength+9]] <- predictmeansBarPlot
		# outputlist[[listlength+10]] <- p_valueMatrix
		# outputlist[[listlength+11]] <- meanTable
        # if (is.null(permlist) || permlist%in%c("NULL", "")) {
          # names(outputlist) <- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                   # "LSD", paste("Pairwise comparison p-value (adjusted by", sQuote(adj), "method)"), 
								   # atvar.levels, "Back Transformed Means", "predictmeansPlot", "predictmeansBKPlot", 
								   # "predictmeansBarPlot", "p_valueMatrix", "mean_table")[1:length(outputlist)]												 
        # }else{
          # names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                   # "Approximated LSD"), paste("Pairwise", sQuote(nsim), "times permuted p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", 
                                                              # paste(sQuote(resvar),collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")), 
															  # "Back Transformed Means with an Approximated 95% CIs (Mean +/- 2*SE)", "predictmeansPlot", "predictmeansBKPlot", 
															  # "predictmeansBarPlot", "p_valueMatrix", "mean_table")[1:length(outputlist)]                            
        # }          
		  # class(outputlist) = "pdmlist"
          # return(outputlist)	        
      # }
    # }else {
      # outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                  # "Standard Error of Differences" = SED.out, LSD = LSD,
                  # "Back Transformed Means" = bkmt, predictmeansPlot=predictmeansPlot, 
				  # predictmeansBKPlot=predictmeansBKPlot, predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable)
	# class(outputlist) = "pdmlist"
          # return(outputlist)	
    # }
  } # else {
    if (pairwise) {
      if (is.null(atvar) || atvar%in%c("NULL", "")) {
        if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj %in% c("none", "bonferroni"))) {
          outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,ndecimal+1), 
					  "Pairwise p-value" = round(t.p.valuem, 4), predictmeansPlot=predictmeansPlot,
					  predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable, p_valueMatrix=p_valueMatrix)
		  class(outputlist) = "pdmlist"
          return(outputlist)			  
        }else{
          outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4), 
					  predictmeansPlot=predictmeansPlot, predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable,
					  p_valueMatrix=p_valueMatrix)
	      class(outputlist) = "pdmlist"
          return(outputlist)
		 }        
      }else{
        outputlist <- vector("list", listlength+8)
        outputlist[[1]] <- mean.table
        outputlist[[2]] <- se.table
        outputlist[[3]] <- SED.out
        outputlist[[4]] <- LSD
		outputlist[[5]] <- paste("For variable", paste(sQuote(resvar), collapse =" and "), "at each level of", paste(sQuote(atvar), collapse =" and "))                   
		for (i in 6: (listlength+5))  outputlist[[i]] <- pmlistTab[[i-5]]
		outputlist[[listlength+6]] <- meanTable
		outputlist[[listlength+7]] <- predictmeansPlot
		outputlist[[listlength+8]] <- predictmeansBarPlot
		outputlist[[listlength+9]] <- p_valueMatrix		
		
		# print(outputlist)
        if (is.null(permlist) || permlist%in%c("NULL", "")) {
          names(outputlist)<- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                   "LSD", paste("Pairwise comparison p-value (adjusted by", sQuote(adj), "method)"), 
								   atvar.levels, "mean_table", "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[1:length(outputlist)]								   
        }else{
          names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                   "Approximated LSD"), paste("Pairwise", sQuote(nsim), "times permuted p-value (adjusted by", sQuote(adj), "method)"), 
								    atvar.levels,								   
								   
								   #paste("Pairwise", sQuote(nsim), "times permuted p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", 
								  # paste(sQuote(resvar), collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")), 
								   "mean_table", 
								   "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[1:length(outputlist)]                     
        }
	  class(outputlist) <- "pdmlist"
      return(outputlist)      
      }
    }else {
      outputlist=list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                  "Standard Error of Differences" = SED.out, LSD = LSD, predictmeansPlot=predictmeansPlot,
				  predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable)
	  class(outputlist) <- "pdmlist"
      return(outputlist)
    }
#  }
}
