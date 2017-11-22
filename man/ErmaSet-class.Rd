\name{ErmaSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ErmaSet-class}
\alias{cellTypes,ErmaSet-method}
\alias{makeErmaSet}
\alias{cellTypes}
\alias{subsetByRanges}

\title{Class \code{"ErmaSet"}}
\description{
Wrap a GenomicFiles instance representing roadmap bed files.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ErmaSet", ...)}.
Simple extension of GenomicFiles.
}
\section{Slots}{
  \describe{
    \item{\code{files}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{rowRanges}:}{Object of class \code{"GenomicRangesORGRangesList"} ~~ }
    \item{\code{colData}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{assays}:}{Object of class \code{"Assays"} ~~ }
    \item{\code{NAMES}:}{Object of class \code{"character_OR_NULL"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\link[GenomicFiles]{GenomicFiles-class}"}, directly.
Class \code{"\link[SummarizedExperiment]{RangedSummarizedExperiment-class}"}, by class "GenomicFiles", distance 2.
Class \code{"\link[SummarizedExperiment]{SummarizedExperiment-class}"}, by class "GenomicFiles", distance 3.
Class \code{"\linkS4class{Vector}"}, by class "GenomicFiles", distance 4.
Class \code{"\linkS4class{Annotated}"}, by class "GenomicFiles", distance 5.
}
\section{Methods}{
  \describe{
    \item{cellTypes}{\code{signature(x = "ErmaSet")}: ... }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}

\usage{
makeErmaSet()
subsetByRanges( ermaset, range )
}
\arguments{
 \item{ermaset}{an instance of \code{\link{ErmaSet-class}}}
 \item{range}{an instance of \code{\link[GenomicRanges]{GRanges-class}}}
}


\note{
We would like to have methods for AnnotationHub entities too.
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
