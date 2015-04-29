\name{genemodel}
\alias{genemodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create GRanges instance with model for a gene
}
\description{
create GRanges instance with model for a gene
}
\usage{
genemodel(sym, genome = "hg19", annoResource = Homo.sapiens, getter = exonsBy, byattr = "gene")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym}{
symbol used as key into Homo.sapiens with keytype SYMBOL
}
  \item{genome}{
reference genome identifier
}
  \item{annoResource}{
OrganismDb instance
}
  \item{getter}{
function that extracts addresses for model from OrganismDb
}
  \item{byattr}{
by element for extraction
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a \code{\link[GenomicRanges]{GRanges}} instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
genemodel("IL33")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
