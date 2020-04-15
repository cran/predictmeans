\name{permmodels}
\alias{permmodels}
\title{
Permutation Test of Linear Model 
}
\description{
This function provides permutation t-tests for coefficients of (fixed) effects and permutation F-tests 
for the terms in a linear model such as \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
\usage{
permmodels(model, data, block=NULL, group=NULL, covariate=NULL, nsim=4999, check=FALSE, 
  exact=FALSE, fo=NULL, prt=TRUE, seed)
}
\arguments{
  \item{model}{
Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. 
}
  \item{data}{
Data object used in the \code{model} fitting.
}
  \item{block}{
Name (in "quotes") for the blocking factor in the \code{data}.
}
  \item{group}{
Name (in "quotes") for the group factor in the \code{data}.
}
  \item{covariate}{
Name(s) (in "quotes") for the covariate variable(s) in the \code{model}.
}
  \item{nsim}{
The number of permutations. The default is 4999.  
}
  \item{check}{
a logical variable to indicate whether or not you want to check permutation results, the default is \code{FALSE}.
}
  \item{exact}{
A logical variable to indicate whether or not exact no. of permutations will be used (applicable only to free the permutation case). The default is FALSE.
}
  \item{fo}{
A model formula used in the \code{model}; \code{fo!=NULL} when the formula is specified by function \code{formula}.}
  \item{prt}{
A logical variable to indicate whether or not to print output on the screen. The default is TRUE.}
  \item{seed}{
Specify a random number generator seed, for reproducible results.
}
}
\value{
The function produces permutation t-test table for coefficients of (fixed) effects, permutation ANOVA table for model terms and a model parameter list \code{permlist}, a list containing \code{nsim=4999} times permutation refitted \code{model} parameters which are used in functions \code{predictmeans} and \code{contrastmeans}. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
## Not run for simplifying process of submiting pkg to CRAN
#library(predictmeans)
#Oats$nitro <- factor(Oats$nitro) 
#fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
## library(lme4)
## fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#
## Permutation Test for model terms
#system.time(
#  permlme <- permmodels(model=fm, data=Oats, block="Block", group="Variety", nsim=999)
#)  
#
## Permutation Test for multiple comparisons
#predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH", 
#  permlist=permlme, plot=FALSE)
#
## Permutation Test for specified contrasts
#cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
}
