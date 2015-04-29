\name{ErmaSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ErmaSet-class}
\alias{makeErmaSet}
\alias{cellTypes,ErmaSet-method}

\title{Class \code{"ErmaSet"}}
\description{
manage epigenomics road map bed file collections, using \code{\link[GenomicFiles]{GenomicFiles-class}} infrastructure
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ErmaSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{files}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{rowData}:}{Object of class \code{"GenomicRangesORGRangesList"} ~~ }
    \item{\code{colData}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{exptData}:}{Object of class \code{"SimpleList"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GenomicFiles}"}, directly.
}
\section{Methods}{
  \describe{
    \item{cellTypes}{\code{signature(x = "ErmaSet")}: return vector of
cell types managed in the object }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
simple extension of GenomicFiles
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("ErmaSet")
makeErmaSet()
}
\keyword{classes}
