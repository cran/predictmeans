\name{Drug}
\alias{Drug}
\docType{data}
\title{\code{Drug} data}

\description{The data is for the comparison of the effectiveness of three analgesic drugs 
to a standard drug, morphine (Finney, Probit analysis, 3rd Edition 1971, p.103). 14 groups 
of mice were tested for response to the drugs at a range of doses. 
}
\usage{
data(Drug)
}
%- maybe also `usage' for other objects documented here.

\format{\code{Drug} is a 14 row data frame with the following columns 
\describe{
\item{Drug}{type of drug.}
\item{Dose}{dose volumn.}
\item{N}{total number of mice in each group.}
\item{R}{number responding mice in each group.}
\item{log10Dose}{log10 transformed dose volumn.}
}
} 

