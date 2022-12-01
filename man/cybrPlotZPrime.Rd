\name{cybrPlotZPrime}
\alias{cybrPlotZPrime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Smoothed Z-scores
%%  ~~function to do ... ~~
}
\description{Plots ggplot geom_line() faceted by chromosome for each factor included of smoothed scores.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrPlotZPrime(zprimedf,
               columns = c("Bulk_Zprime", "Parent_Zprime", "Interaction_Zprime"),
               chromosomes = "All",
               title = "Smoothed Z Scores",
               yeast = TRUE,
               colvalues = c("#345F6F", "#D7335C", "#FFB05C")
               )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zprimedf}{Dataframe including data to plot; should be smoothed Z-scores, but will plot a line plot for any data. Recommended output of cybrSmoothBSAWindows(), but can take in cybrBSA_GLM.
%%     ~~Describe \code{x} here~~
}
  \item{columns}{Character vector of column which should be plotted. Recommended to plot smoothed scores as lines, rather than multiple unsmoothed
%%     ~~Describe \code{x} here~~
}
  \item{chromosomes}{Default "All" will include all chromosomes besides "M" (mitochondria); for any other subset, list as character vector
%%     ~~Describe \code{x} here~~
}
  \item{title}{Character string title for plot
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{Boolean; if TRUE will convert chromosomes to leveled factor; if FALSE will leave original chromosome names
%%     ~~Describe \code{x} here~~
}
  \item{colvalues}{Character vector of colors to plot
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot of each chromosome and the Z scores of each position connected by lines
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Y Buzby
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
