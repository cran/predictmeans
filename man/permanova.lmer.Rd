\name{permanova.lmer}
\alias{permanova.lmer}
\title{
Permutation ANOVA for \code{lmer} Model
}
\description{
This function provides permutation ANOVA for \code{lmer} model.
}
\usage{
permanova.lmer(model, nperm = 999, ncore=3L, type = c("I", "II", "III",  "1", "2", "3"),
 ...) 
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
  \item{type}{
The type of ANOVA table requested (using SAS terminology) with Type I being the familiar 
sequential ANOVA table.
}
  \item{...}{
Use to setup option: seed -- Specify a random number generator seed, for reproducible 
results.
}
}
\value{
Permutation ANOVA table. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
## NOT RUN
# library(predictmeans)
# Oats$nitro <- factor(Oats$nitro) 
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# 
# # Permutation Test for model terms
# permanova.lmer(fm)
# permanova.lmer(fm, type=2)
# 
# # Compare to F test
# fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# anova(fm0)
}
