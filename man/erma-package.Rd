\name{erma-package}
\alias{erma-package}
\alias{erma}
\docType{package}
\title{
\packageTitle{erma}
}
\description{
\packageDescription{erma}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{erma}
\packageIndices{erma}

This package provides infrastructure for working with products
of the NIH epigenome roadmap project.
}
\author{
\packageAuthor{erma}

Maintainer: \packageMaintainer{erma}
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
data(farhSE)
farhSE
ee = makeErmaSet()
stateProfile(ee[,1:4], symbol="A2M", orient5p=TRUE)
}
