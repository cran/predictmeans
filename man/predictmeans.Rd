\name{predictmeans}
\alias{predictmeans}
\title{
Predicted Means of a Linear Model
}
\description{
  This function obtains predicted means, SE of means, SED of means, LSDs and plots of means
  with SE bar or LSD bar for parametric models such as \code{aov}, \code{lm}, 
  \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. The function also perfomrs pairwise comparisons
  and permutation tests.
}
\usage{
predictmeans(model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL,
  lsd_bar=TRUE, level=0.05, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I,
  transOff = 0, responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, plottitle=NULL, 
  plotxlab=NULL, plotylab=NULL, mplot=TRUE, barplot=FALSE, pplot=TRUE, bkplot=TRUE, 
  plot=TRUE, jitterv=0.2, basesz=12L, prtnum=TRUE, prtplt=TRUE, newwd=FALSE, 
  permlist=NULL, ncore=3L, ndecimal=4L)
}
\arguments{
  \item{model}{
  Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
  \item{modelterm}{
  Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
  The \code{modelterm} must be factors and given exactly as it appears in the printed model, e.g. "A" or "A:B". 
}
  \item{data}{
  In some cases, you need to provide the data set used in model fitting, especially when you have applied some 
  variable trnasformation in the model.}
  \item{pairwise}{
  An option for showing pair-wise LSDs and p-values, or not. The default is FALSE. 
}
  \item{atvar}{
  When \code{pairwise = TRUE}, a quoted name indicating within levels of which variable in 
  \code{modelterm} the multiple comparison will be performed.
}
  \item{adj}{
  Name (in "quote") for indicating a method for adjusting p-values of pairwise comparisons. 
  The choices are "none", "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr". 
  The default method is "none". Note that LSD can't be adjusted except for "bonferroni" method.
}
  \item{Df}{
  A degree of freedom for calculating LSD. For the above models, Df is obtained from the function automatically.
}
  \item{lsd_bar}{
  A logical variable to indicate to print an average LSD or SED bar on the means plot. The default is TRUE.
}
  \item{level}{
  A significant level for calculating LSD, CI etc. The default value is 0.05.
}
  \item{covariate}{
  A numerical vector to specify values of covariates for calculating predicted means. The default 
  values are the means of the associated covariates.
}
  \item{meandecr}{
  A logical variable to indicate whether to print letters for multiple comparisons by decreasing order of means in the mean_table. 
  The default is NULL which indicates the mean order follows the associated factor levels.
}
  \item{letterCI}{
  A logical variable to indicate printed letters for multiple comparisons by whether or not CI overlap in the mean_table. 
  The default is FALSE. Note that the method of examining overlap is more conservative (i.e., rejects the null hypothesis less often) than the standard method when the null hypothesis is true.
}
  \item{trans}{
  A function object for calculating the back transformed means, e.g. \code{trans=exp}.
}
  \item{transOff}{
When you use \code{trans=exp(x+1)}, then \code{transOff=1}, the default is 0.
}
  \item{responsen}{
  Name (in "quotes") of the back transformed response variable in the \code{model}.
}
  \item{count}{
  An option for indicating the back transformed mean values are counts or not. The default is FALSE.
} 
 \item{plotord}{
  A numeric (or character) vector specifying the order of plotting for two or three way interaction (e.g. 
  \code{plotord = c(2, 1, 3)} or \code{plotord = c("B", "A", "C")} will put the second variable 'B' in 
  \code{modelterm} on the \code{X} axis, the first variable 'A' as the grouping variable, and the third one 
  'C' as the panel variable). The defaults are \code{c(1, 2)} and \code{c(1, 2, 3)} for two and three way 
  interactions. 
} 
  \item{lineplot}{
  An option for drawing a line chart, or dot chart. The default is TRUE.
}
  \item{plottitle}{
  A character vector specifying the main title for plot(s). The default is NULL. 
}
  \item{plotxlab}{
  A character vector specifying the x label for plot(s). The default is NULL. 
}
  \item{plotylab}{
  A character vector specifying the y label for plot(s). The default is NULL. 
}
  \item{mplot}{
  An option for drawing a means plot, or not. The default is TRUE.
} 
  \item{barplot}{
  An option for drawing a bar chart, or not. The default is FALSE.
}
  \item{pplot}{
  An option for drawing a p-values plot, or not when there are more than six p-values. The default is TRUE.
}
  \item{bkplot}{
  An option for drawing back transformed plot, or not. The default is TRUE.
}
  \item{plot}{
  An option for drawing plots, or not. The default is TRUE.
}
  \item{jitterv}{
  A degree of jitter in x and y direction in the back transformed means graph. The default is zero.
}
  \item{basesz}{
  The base font size. The default is 12.
}
  \item{prtnum}{
  An option for printing covariate information on the screen, or not. The default is TRUE.
}
  \item{prtplt}{
  An option for printing plots on the screen, or not. The default is TRUE.
}
  \item{newwd}{
  A logical variable to indicate whether to print graph in a new window. The default is FALSE.
}
  \item{permlist}{
  A model parameter list produced by the function \code{permmodels}. When \code{permlist != NULL}, 
  the option \code{Df} will be non-functional. This is a key option for pairwise comparisons via 
  permutation tests.
} 
  \item{ncore}{
Number of core for parallel computing when \code{permlist != NULL}, the default value is 3.
} 
  \item{ndecimal}{
  An option for specifying number of decimal point to be print at predicted means table. The default is 4.
}
}

\value{
 \item{Predicted Means}{A table of predicted means.}
 \item{Standard Error of Means}{A table of standard errors of predicted means.}
 \item{Standard Error of Differences}{Standard errors of differences between predicted means.}
 \item{LSD}{Least significant differences between predicted means.}
 \item{Pairwise p-value}{A matrix with t-values above the diagonal and p-values below the diagonal, or 
 matrix of pairwise comparison p-values for each level of \code{atvar}.}
 \item{mean_table}{A summary of predicted means result including 'Predicted means', 'Standard error', 'Df' 
 and 'CIs'. When \code{trans!=NULL} or \code{trans!=I}, a table of back transformed means with CIs are also shown.}
 \item{predictmeansPlot}{ggplot of predicted means.}
 \item{predictmeansBKPlot}{ggplot of back transformed means.}
 \item{predictmeansBarPlot}{gg bar plot of predicted means.}
 \item{p_valueMatrix}{p_value matrix for pairwise comparison.}
}

\note{
  The \code{predictmeans} function becomes confused if a factor or covariate is changed to the other
  in a model formula. Consequently, formulae that include calls \code{as.factor}, \code{factor}, or \code{numeric}
  (e.g. \code{as.factor(income)}) will cause errors. Instead, create the modified variables outside of the model
  formula (e.g., \code{fincome <- as.factor(income)}) and then use them in the model formula.
  
  Factors cannot have colons in level names (e.g., \code{"level:A"}); the \code{predictmeans} function will confuse the
  colons with interactions; rename levels to avoid colons.
  
  For \code{predictmeans} function, it is assumed that methods \code{coef}, \code{vcov}, \code{model.matrix}, \code{model.frame} and \code{terms} are available for \code{model}. 
}

\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\references{
	Maghsoodloo Saeed, Ching-Ying Huang (2010),
	\emph{Comparing the overlapping of two independent confidence intervals with a single confidence interval for two normal population parameters}, Journal of Statistical Planning and Inference, 140(11), 3295-3305. https://www.sciencedirect.com/science/article/pii/S0378375810002405. 
	
	Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
    \emph{Simultaneous Inference in General Parametric Models. Biometrical}, Journal 50(3), 346-363.
  
	Welham S., Cullis B., Gogel B., Gilmour A., & Thompson R. (2004),
	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347.
}

\examples{
  library(predictmeans)
  ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
  Oats$nitro <- factor(Oats$nitro)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  predictmeans(fm, "nitro", adj="BH")
  predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", line=FALSE)
  predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", 
    barplot=TRUE, line=FALSE)
  names(predictout)
  print(predictout$predictmeansPlot)
  print(predictout$predictmeansBarPlot)
}