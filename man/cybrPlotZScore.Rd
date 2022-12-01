\name{cybrPlotZScore}
\alias{cybrPlotZScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Z-scores
%%  ~~function to do ... ~~
}
\description{Uses ggplot geom_point() to plot Z scores faceted by chromosome. Requires ggplot2 and dplyr.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrPlotZScore(zprimedf,
               column = "Bulk_Z",
               chromosomes = "All",
               title = "Z Scores by Position",
               yeast = TRUE,
               dotcolor = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zprimedf}{Dataframe of Z scores; can be output of cybrSmoothBSAWindows() or cyberBSA_GLM()
%%     ~~Describe \code{x} here~~
}
  \item{column}{Character vector of column which should be plotted
%%     ~~Describe \code{x} here~~
}
  \item{chromosomes}{Default "All" will include all chromosomes besides "M" (mitochondria); for any other subset, list as character vector
%%     ~~Describe \code{x} here~~
}
  \item{title}{Character value for title of plot
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{Boolean, if TRUE will convert chromosomes to leveled factor, if FALSE will not change
%%     ~~Describe \code{x} here~~
}
  \item{dotcolor}{Character value for a single color for entire plot
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cybrPlotZPrime(CuSO4_wholegenomeBSA, columns = c("Bulk_Zprime"))

cybrPlotZScore(CuSO4_wholegenomeBSA)

cybrPlotZPrime(wholegenomeBSA, chromosomes = ChromKey$chromosomes[2:16])

cybrPlotZScore(wholegenomeBSA)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
