\name{semipred}
\alias{semipred}
\title{
  Predicted Means of a Semi Paramatric Model with Covariate Variable(s)
}
\description{
  This function produces predicted means with graph for a semi paramatric model with new set of covariate values.
}
\usage{
  semipred(semireg, modelterm=NULL, covariate, sm_term=NULL, contr=NULL,
           covariateV=NULL, boundary=NULL, level=0.05, trans=NULL, trellis=TRUE, 
		   scales=c("fixed", "free", "free_x", "free_y"),
		   plotord=NULL, ci=TRUE, point=TRUE, jitterv=0, threeD=FALSE, prt=TRUE)
}
\arguments{
  \item{semireg}{
  A list object returned by \code{semireg}. 
}
  \item{modelterm}{
  Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
  The \code{modelterm} must be given exactly as it appears in \code{semireg} model, e.g. "A" or "A:B". In case \code{modelterm} is the same as \code{covariate} or \code{NULL}, them \code{semipred} will product predictmeans with CI based on \code{covariate} only with out any grouping.
}
  \item{covariate}{
  Name (in "quotes") of one or two (for \code{Ztps} smooth) the covariate variables in the \code{semireg} shuch as "x1" or c("x1", "x2").
}
  \item{sm_term}{
Names (in "quotes") of smooth terms (from \code{smoothZ} list in \code{semireg} model) used in the prediction such as "sm1_grp" or c("sm1_grp", "sm2_grp"). The default is using all smooth terms which is \code{sm_term=NULL}. 
} 
  \item{contr}{
A numeric (or character) vector with length of two (e.g. c(4, 1) or c("d", "a")) which indicates to produce predicted mean with CI for difference between \code{modelterm} level 4 ("d") vs level 1 ("a") along \code{covariate}.
} 
  \item{covariateV}{
  A numeric vector or matrix, then semipred will produce the result for \code{covariate} at value of \code{covariateV}.
}
  \item{boundary}{
  A matrix or data frame of two columns, used to specify boundary of longitude and latitude, it is functional when the length of covariate is two.
}
  \item{level}{
  A significant level for calculating confident interval. The default value is 0.05.
}
  \item{trans}{
A function object for calculating the back transformed means, e.g. \code{trans=exp}.
}  
  \item{trellis}{
  A logical scalar. If set to TRUE (default), a trellis plots of predicted means with CI will be drawn.
}
  \item{scales}{
  Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y") in a trellis graph?
}
  \item{plotord}{
  A numeric vector specifying the order of plotting for two way interaction (e.g. 
  \code{plotord = c(2, 1)} will put the second variable in \code{modelterm} on the \code{X} 
  axis, the first variable as the grouping variable, and the third one as the panel variable). 
  The defaults are \code{c(1, 2)} for two way interactions. 
}
  \item{ci}{
A logical scalar to indicate whether to print confidence interval. The default value is TRUE.
} 
\item{point}{
A logical scalar to indicate whether to print raw data points. The default value is TRUE.
} 
\item{jitterv}{
  A degree of jitter in x and y direction in the graph. The default is zero.
}
  \item{threeD}{
A logical scalar to indicate whether to produce a 3-D plot or not. The default value is FALSE.
} 
  \item{prt}{
A logical scalar to indicate whether to produce plots on the screen or not. The default value is TRUE.
}
}

\value{
\item{plt}{A ggplot object.}
\item{pred_df}{A data.frame with predcted data.}
}

\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
## NOT RUN 
# library(predictmeans)
# data(Dialyzer, package="nlme")
# help(Dialyzer)
# str(Dialyzer)
# 
# library(ggplot2)
# ggplot(Dialyzer, aes(x=pressure, y=rate, col=QB)) +
#   geom_line() +
#   facet_wrap(vars(Subject))
# 
# fm <- semireg(rate ~ pressure*QB+(pressure|Subject),
#               smoothZ=list(
#                 qb_grp=smZ(pressure, k=4, by=QB, group=TRUE)
#               ),
#               data=Dialyzer)
# str(fm$data)
# summary(fm$semer)
# residplot(fm$semer, group="QB")
# anova(fm$semer)
# ranova(fm$semer)
# R2_glmm(fm$semer)
# ap_out1 <- semipred(fm, "QB", "pressure")
# str(ap_out1$pred_df)
# ap_out2 <- semipred(fm, "QB", "pressure", contr=c(1,2))
# str(ap_out2$pred_df)
# 
# data(sleepstudy, package="lme4")
# help(sleepstudy)
# str(sleepstudy)
# library(latticeExtra)
# xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
#        layout = c(9, 2), type = c("g", "p", "r"),
#        index.cond = function(x, y) coef(lm(y ~ x))[2],
#        xlab = "Days of sleep deprivation",
#        ylab = "Average reaction time (ms)",
#        as.table = TRUE)
# 
# sleep.semi <- semireg(Reaction ~ Days*Subject,
#                       smoothZ=list(
#                         sub_grp=smZ(Days, by=Subject, group=TRUE)
#                       ),
#                       data=sleepstudy)
# residplot(sleep.semi$semer)
# summary(sleep.semi$semer)
# anova(sleep.semi$semer)
# ranova(sleep.semi$semer)
# R2_glmm(sleep.semi$semer)
# 
# predout1 <- semipred(sleep.semi, "Subject", "Days")
# str(predout1$pred_df)
# predout2 <- semipred(sleep.semi, "Subject", "Days", contr = c(6,1))
# str(predout2$pred_df)
}
