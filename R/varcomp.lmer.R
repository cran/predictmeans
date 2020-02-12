varcomp.lmer <- function(model, level=0.95) {
  if (!inherits(model, "merMod")) stop("The model must be a lmer object!")
  lmemod.ML <- refitML(model)
  varcomp <- data.frame(VarCorr(model))
  
  Rtms <- gsub("_\\(Intercept\\)", "", gsub("_NA", "",  do.call("paste", c(varcomp[, !names(varcomp)%in%c("vcov", "sdcor"), drop=FALSE], sep="_"))))
  
  lmemod.MLvcov <- waldVar2(lmemod.ML)
  dimnames(lmemod.MLvcov) <- list(Rtms,Rtms)
  varcomp$Std.error <- 2*sqrt(varcomp$vcov)*sqrt(diag(lmemod.MLvcov))
  pars <- as.data.frame(VarCorr(model),order="lower.tri")[,"sdcor"]
  names(pars) <- Rtms
  vhack2 <- list(coefficients=pars, vcov=lmemod.MLvcov)
  wci <- confintlmer(vhack2, level=level)
  
  wci[wci<0] <- 0
  varcomp <- round(cbind(varcomp[, c("vcov", "Std.error")], wci^2), 4)
  return(varcomp)
}
