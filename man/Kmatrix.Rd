\name{Kmatrix}
\alias{Kmatrix}
\title{
Matrix of Coefficients in a Linear Model
}
\description{
  This function obtains a matrix of coefficients for parametric models such as \code{aov}, \code{lm}, 
  \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
\usage{
Kmatrix(model, modelterm, covariate=NULL, as.is=FALSE, prtnum=FALSE)
}
\arguments{
  \item{model}{
Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
  \item{modelterm}{
  Name (in "quotes") for indicating which model term's predicted mean to be calculated. 
  The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B". 
}
  \item{covariate}{
A numerical vector to specify values of covariates for calculating predicted means, 
default values are the means of the associated covariates. It also can be the name of covariate in the model.
}
  \item{as.is}{
A logic value to specify wheather or not using original covariate values, the default is FALSE.
}
  \item{prtnum}{
An option for printing covariate info on the screen or not. The default is FALSE.
}
}

\value{
 \item{K}{Coefficients matrix}
 \item{fctnames}{A model frame contains factor(s) info in the model.}
 \item{response}{The name of response variable in the model.}
}
\references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347. 
}
\author{
This function heavily depends on the codes from package "lsmeans".
}

\examples{
  library(predictmeans)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  Kmatrix(fm, "Variety", prtnum=TRUE)
  Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)
#  Kmatrix(fm, "Variety", "nitro")
}

