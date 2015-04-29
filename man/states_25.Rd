\name{states_25}
\alias{states_25}
\docType{data}
\title{
metadata on states of 25-state model of chromatin from ChromImpute
}
\description{
metadata on states of 25-state model of chromatin from ChromImpute
}
\usage{data("states_25")}
\format{
  A data frame with 25 observations on the following 5 variables.
  \describe{
    \item{\code{STATENO.}}{a numeric vector}
    \item{\code{MNEMONIC}}{a character vector}
    \item{\code{DESCRIPTION}}{a character vector}
    \item{\code{COLOR.NAME}}{a character vector, partly non-compliant with R colors}
    \item{\code{COLOR.CODE}}{a character vector, RGB numerics, comma-delimited}
    \item{\code{rgb}}{a character vector, RGB scores in R atomic format}
  }
}
\details{
Some modifications needed to remove registered trademark symbol
mistakenly present in 'regulatory' and
non-ascii elements of prime notations
}
\source{
retrieved from \url{http://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp} 28 April 2015
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(states_25)
## maybe str(states_25) ; plot(states_25) ...
}
\keyword{datasets}
