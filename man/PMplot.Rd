\name{PMplot}
\alias{PMplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Level Plot of a Matrix of p-values.
}
\description{
Creates a plot of p-values of pairwise comparisons.
}
\usage{
PMplot(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73, 
  newwd=TRUE)  
}

\arguments{
  \item{pmatrix}{
  A matrix with p-values from pairwise comparisons. (This is a lower triangle matrix.)
}
  \item{level}{
  The level of p-value to be highlighted. Default is 0.05.
}
  \item{mtitle}{
  The main title in the graph.
}
  \item{xylabel}{
  The x and y labels in the graph. 
}
  \item{margin}{
  A value for specifying x and y margins in the graph. The default value is 5.
}
  \item{legendx}{
  A value for specifying x coordinate of legend. The default value is 0.73.
}
\item{newwd}{
A logical variable to indicate whether to print graph in a new window. The default is TRUE.
}
}

\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
library(predictmeans)
set.seed(2013)
pvalues <- runif(28)
pmatrix <- matrix(0,8,8)
pmatrix[lower.tri(pmatrix)] <- pvalues
round(pmatrix, 4)
PMplot(pmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
