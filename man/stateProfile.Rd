\name{stateProfile}
\alias{stateProfile}
\alias{stateProf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a ggplot2 visualization of states over genomic addresses
for a family of cell types managed in an ErmaSet instance
}
\description{
Create a ggplot2 visualization of states over genomic addresses
for a family of cell types managed in an ErmaSet instance.
}
\usage{
stateProfile(ermaset, symbol = "IL33", width = 50000, ctsize = 7, shortCellType = TRUE)
stateProf(ermaset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ermaset}{
instance of \code{\link[erma]{ErmaSet-class}}
}
  \item{symbol}{
gene symbol resolvable in \code{\link[Homo.sapiens]{Homo.sapiens}}
}
  \item{width}{
width of region upstream of gene to be visualized
}
  \item{ctsize}{
font size for cell type labels
}
  \item{shortCellType}{
logical, if FALSE, full cell type labels are used, otherwise
ad hoc abbreviations are used for facet labeling
}
}
\details{
stateProf is a demonstration shiny app calling stateProfile
}
\value{
instance of c("gg", "ggplot")
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
}

%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
ermaset = makeErmaSet()
stateProfile(ermaset, "CD28")
}
\keyword{ models }
