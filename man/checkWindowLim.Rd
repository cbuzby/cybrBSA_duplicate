\name{checkWindowLim}
\alias{checkWindowLim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check maximum window per chromosome
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculate maximum window size W to use in cybrSmoothBSAWindows(window = W) per chromosome; gives the number of SNPs or datapoints and the minimum window size that can be used for smoothing.
}
\usage{
checkWindowLim(Dataset, includechr = TRUE, exceptchr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dataset}{
%%     ~~Describe \code{Dataset} here~~
    Data frame with CHROM and POS columns
}
  \item{includechr}{
%%     ~~Describe \code{includechr} here~~
    If TRUE, returns all chromosomes; if FALSE or if exceptchr is not NULL, returns the CHROM column that matches the includechr vector, and excludes those in the exceptchr vector.
}
  \item{exceptchr}{
%%     ~~Describe \code{exceptchr} here~~
    If is.null = FALSE, excludes those in vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
Dataframe containing three columns:
%%  If it is a LIST, use
 \item{CHROM }{Unique chromosomes in CHROM input}
 \item{SNPs }{Number of unique SNPs per chromosome}
 \item{maxW }{Half the number of unique SNPs per chromosome; this is the maximum window size for this chromosome}

%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Cassandra Buzby
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
This function should be used to determine which chromosomes should be excluded when running cybrSmoothBSAWindows() to prevent errors. Updated versions will return smoothed values for those with the min number of datapoints, but that is yet to be implemented.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

sampledata <- data.frame(CHROM = factor(sample(ChromKey$chromosomes,
                                        size = 100,
                                        replace = TRUE,
                                        prob = c(1,3,1.2,10,5,5,5,5,4,
                                                 4,3,6,4,3,2,5,0.5))),
                          POS = sample(1:5000, 1000, replace = TRUE),
                          DummyData = sample(1:5000, 1000, replace = TRUE))

checkWindowLim(sampledata)

checkWindowLim(sampledata, includechr = c("I", "II", "III", "IV", "V"))

#Find the maximum window size for the entire dataset besides M and Chr III
min(checkWindowLim(sampledata, exceptchr = c("M", "III"))$maxW)

#Find which chromosomes would not work with a window size of 20
checkWindowLim(sampledata)$CHROM[which(checkWindowLim(sampledata)$maxW < 20)]

-------------------------------------------------------------------------------
# The function is currently defined as:
-------------------------------------------------------------------------------

function (Dataset, includechr = TRUE, exceptchr = NULL)
{
    if (is.null(exceptchr) == FALSE) {
        includechr = FALSE
    }
    if (includechr != TRUE) {
        if (is.null(exceptchr)) {
            Dataset <- Dataset \%>\% filter(CHROM \%in\% includechr)
        }
        else {
            Dataset <- Dataset \%>\% filter(CHROM \%in\% exceptchr ==
                FALSE | CHROM \%in\% includechr)
        }
    }
    tablecounts <- Dataset \%>\% group_by(CHROM) \%>\% summarise(CHROM = CHROM,
        SNPs = length(unique(POS))) \%>\% distinct() \%>\% mutate(minW = SNPs/2)
    return(tablecounts)
  }
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
