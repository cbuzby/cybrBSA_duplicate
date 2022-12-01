\name{cybrQualityFilter}
\alias{cybrQualityFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to filter for quality and remove positions that do not contain all values
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrQualityFilter(gatkdf, GQcutoff = 98, cleandata = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gatkdf}{Input dataframe, takes direct output from cybrInputGATKTable()
%%     ~~Describe \code{x} here~~
}
  \item{GQcutoff}{Numeric input between 0 and 100; 98 recommended for standard GATK output
%%     ~~Describe \code{x} here~~
}
  \item{cleandata}{Boolean, if TRUE will remove positions that do not contain all categories to prevent downstream errors in logistic regression; if FALSE will return only the values where Genome Quality is higher than GQcutoff input value
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Dataframe }{Includes all input columns but only rows which pass filters}
%%  \item{comp2 }{Description of 'comp2'}
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
mydatatotest = "HGV.SortedCat.vcf.output.table"

mydf <- cybrInputGATKTable(mydatatotest)

qualitydf <- cybrQualityFilter(mydf)

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

testmerge <- cybrIDAlleles(BSAdfstart = qualitydf, Parentdf = parentSNPids, yeast = TRUE) %>% na.omit()
}

