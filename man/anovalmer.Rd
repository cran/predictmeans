\name{anovalmer}
\alias{anovalmer}
\title{
ANOVA of a Linear Mixed Effects Model produced by \code{lmer} function
}
\description{
This function provides ANOVA of a linear mixed effects model for \code{lmer} object with model terms' denominate degree of freedom estimated by Kenward-Roger method.
}
\usage{
anovalmer(model, DDf=NULL)
}
\arguments{
  \item{model}{
Model object returned by \code{lmer}. 
}
  \item{DDf}{
Model terms' denominate degree of freedom specified by user.
}
}

\value{
ANOVA table. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
## Not run for simplifying process of submiting pkg to CRAN
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
anovalmer(fm)
}
