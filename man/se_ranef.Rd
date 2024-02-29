\name{se_ranef}
\alias{se_ranef}

\title{Extract Standard Errors of Model Random Effects}

\description{
   These functions extract standard errors of model random effects
   from objects returned by modeling functions.
}
\usage{
  se_ranef(object, rand_term=NULL)
}

\arguments{
  \item{object}{object of \code{merMod} and \code{glmmTMB} fit}
  \item{rand_term}{a name of random term in the model.}
}
\value{
  \code{se_ranef} gives a list of standard errors for \code{ranef}.
}
\details{
  \code{se_ranef} extracts standard errors of the random effects
    from objects returned by lmer, glmer and glmmTMB functions.
}

\references{This function is modified from function 'se.ranef' at package 'arm'.}

\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}
