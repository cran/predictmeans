\name{permanova.lmer}
\alias{permanova.lmer}
\title{
Permutation ANOVA for \code{lmer} Model
}
\description{
This function provides permutation ANOVA for \code{lmer} model.
}
\usage{
permanova.lmer(model, nperm = 999, ncore=3, drop=TRUE, ...) 
}
\arguments{
  \item{model}{
Model object returned by \code{lmer}. 
}
  \item{nperm}{
Number of permutation, the default value is 999.
}
  \item{ncore}{
Number of core for parallel computing, the default value is 3.
}
  \item{drop}{
Drop or add model term sequentially, the default value is TRUE.
}
  \item{...}{
Use to setup option: seed -- Specify a random number generator seed, for reproducible results.
}
}
\value{
Permutation ANOVA table. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Permutation Test for model terms
# permanova.lmer(fm)
# permanova.lmer(fm, drop=FALSE)
## Compare to F test
# fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# anova(fm0)
}
