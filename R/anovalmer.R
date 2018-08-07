anovalmer <- function(model, DDf=NULL)  {
  if (!inherits(model, "merMod")) stop("The model must be a lmer object!")
  aTable <- anova(model, ddf="Kenward-Roger", type=1)
  if (is.null(DDf) || DDf%in%c("NULL", "")){
  DDf <- p.value <- numeric(nrow(aTable))
  termlabel <- row.names(aTable)
  names(DDf) <- names(p.value) <- termlabel
  
  for (vars in termlabel) {
    KRout <- tryCatch({calcKRDDF(model, vars)}, error = function(e) { NULL })
    if(!is.null(KRout)){
      DDf[vars] <- getKR(KRout, "ddf")
      p.value[vars] <- getKR(KRout, "p.value")
    }else DDf[vars] <- p.value[vars] <- NA    
  }
  if(all(!is.na(DDf))){
    aTable$DDf <- DDf
    aTable$p.value <- round(p.value, 5)
    attr(aTable, "heading") <- "Analysis of Variance Table of type I with Kenward-Roger\napproximation for degrees of freedom"
  }else message("anova from lme4 is returned some error for computing\nKenward-Roger approximation for degrees of freedom")
  return(aTable)
  }else{
    if (length(DDf)!=nrow(aTable)) stop("Please provide suitable DenDF!")
	aTable$DDf <- DDf
	p.value <- 1- pf(aTable[, "F value"], aTable[, "Df"], DDf)
    aTable$p.value <- round(p.value, 5)
	return(aTable)
  }
}

