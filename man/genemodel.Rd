\name{genemodel}
\alias{genemodel}
\alias{geneTxRange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create GRanges instance with model for a gene
}
\description{
create GRanges instance with model for a gene
}
\usage{
genemodel(key, keytype, annoResource = Homo.sapiens)
geneTxRange(sym, annoResource = Homo.sapiens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym}{
symbol used as key into \code{annoResource} with keytype SYMBOL
}
  \item{key}{
string used as key into \code{annoResource} with keytype \code{keytype}
}
  \item{keytype}{
string used as keytype for \code{select} with \code{annoResource} 
}
  \item{annoResource}{
OrganismDb instance
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
\note{
genemodel revised Aug 10 2015.
Direct operations on Homo.sapiens, much
faster.
geneTxRange added Aug 10 2015.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
genemodel("IL33")
geneTxRange("IL33")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
