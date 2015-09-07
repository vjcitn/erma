\name{stateProfile}
\alias{stateProfile}
\alias{stateProf}
\alias{csProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a ggplot2 visualization of chromatin states over genomic addresses
for a family of cell types managed in an ErmaSet instance
}
\description{
Create a ggplot2 visualization of chromatin states over genomic addresses
for a family of cell types managed in an ErmaSet instance.
}
\usage{
stateProfile(ermaset, symbol = "IL33", upstream = 2000, 
   downstream = 200, ctsize = 10, 
   shortCellType = TRUE, tsswidth=3 )
stateProf(ermaset, shortCellType=TRUE, ctsize=10)
csProfile(ermaset, symbol, upstream = 2000, downstream = 200,
                 useShiny = FALSE, ctsize = 10, shortCellType = TRUE,
                 tsswidth = 3)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ermaset}{
instance of \code{\link[erma]{ErmaSet-class}}
}
  \item{symbol}{
gene symbol resolvable in \code{\link[Homo.sapiens]{Homo.sapiens}}
}
  \item{upstream, downstream}{
parameters passed to \code{\link[GenomicRanges]{promoters}} to
limit region to view
}
  \item{ctsize}{
font size for cell type labels
}
  \item{shortCellType}{
logical, if FALSE, full cell type labels are used, otherwise
ad hoc abbreviations are used for facet labeling
}
  \item{useShiny}{logical; if TRUE, shiny app is run permitting
interactive selection of gene, and scope of view}
  \item{tsswidth}{width in base pairs of the base of a black
rectangle used to depict location of transcription start site}
}
%\details{
%stateProf is a demonstration shiny app calling stateProfile
%}
\value{
if useShiny is FALSE, an instance of c("gg", "ggplot") is returned
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
# set useShiny=TRUE for interactive display
csProfile(ermaset[,1:5], "CD28")
}
\keyword{ models }
