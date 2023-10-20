permanova.lmer <- function(model, nperm = 999, ncore=3, type = c("I", "II", "III",  "1", "2", "3"), ...)  { # perms

  if (!inherits(model, "lmerMod")) stop("The model must be a lmer object!")
  if (inherits(try(refit(model, getME(model, "y")), TRUE), "try-error")) stop("Please remove missing values in your model!")
  type <- type[1L]
  if(!is.character(type)) type <- as.character(type)
  type <- match.arg(type)
  if(type %in% c("I", "II"))
    type <- as.character(as.integer(as.roman(type)))
	
  aTable <- anova(model, type=type)
  Perm.p <- numeric(nrow(aTable))
  termlabel1 <- termlabel2 <- row.names(aTable)
  names(Perm.p) <- termlabel1
  
  termlabel0 <- attr(terms(model), "term.labels")
  model.0 <- update(model, as.formula(paste(".~. -", paste(termlabel0, collapse="-"))))

  for (vars in termlabel1) {
    if (type=="3") {
    varsn <- unlist(strsplit(vars, "\\:"))
    for (i in varsn) termlabel2 <- termlabel2[grep(i, termlabel2)]
    termlabel <- paste(termlabel2, collapse="-")
    model.b <- update( model, as.formula(paste(".~. -", termlabel)))
	Perm.p[vars] <- permlmer(model.b, model, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
    termlabel2 <- termlabel1
	}else if (type=="2"){
	varsn <- unlist(strsplit(vars, "\\:"))	
	model.1 <- update(model.0, as.formula(paste(".~.+", paste(termlabel0[attr(terms(model), "order")%in%(1:length(varsn))], collapse="+"))))
	model.2 <- update(model.1, as.formula(paste(".~. -", vars)))
	Perm.p[vars] <- permlmer(model.2, model.1, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]	
	}else if (type=="1"){
	model.1 <- update(model.0, as.formula(paste(".~. +", vars)))
	Perm.p[vars] <- permlmer(model.0, model.1, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
	model.0 <- model.1	
	}else NULL
  }
  aTable$Perm.p <- Perm.p
  return(aTable)
} 