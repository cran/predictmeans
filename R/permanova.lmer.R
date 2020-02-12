permanova.lmer <- function(model, nperm = 999, ncore=3, drop=TRUE, ...)  { # perms
  if (!inherits(model, "merMod")) stop("The model must be a lmer object!")
  aTable <- anovalmer(model)
  Perm.p <- numeric(nrow(aTable))
  termlabel1 <- termlabel2 <- row.names(aTable)
  names(Perm.p) <- termlabel1
  
  termlabel0 <- attr(terms(model), "term.labels")
  model.0 <- update(model, as.formula(paste(".~. -", paste(termlabel0, collapse="-"))))

  for (vars in termlabel1) {
    if (drop) {
    varsn <- unlist(strsplit(vars, "\\:"))
    for (i in varsn) termlabel2 <- termlabel2[grep(i, termlabel2)]
    termlabel <- paste(termlabel2, collapse="-")
    model.b <- update( model, as.formula(paste(".~. -", termlabel)))
	Perm.p[vars] <- permlmer(model.b, model, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
    termlabel2 <- termlabel1
	}else{
	model.1 <- update(model.0, as.formula(paste(".~. +", vars)))
	Perm.p[vars] <- permlmer(model.0, model.1, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
	model.0 <- model.1	
	}
  }
  aTable$Perm.p <- Perm.p
  return(aTable)
}