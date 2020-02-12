\name{varcomp.lmer}
\alias{varcomp.lmer}
\title{
Calculate stder and CI of variance components for \code{lmer} Model
}
\description{
This function calculates stder and CI of variance components for \code{lmer} Model.
}
\usage{
varcomp.lmer(model, level=0.95) 
}
\arguments{
  \item{model}{
Model object returned by \code{lmer}. 
}
  \item{level}{
level of confidence of CI, the default value is 0.95.
}
}
\value{
Variance components table. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
varcomp.lmer(fm)
# Compare with simulation
confint(fm)^2
}
