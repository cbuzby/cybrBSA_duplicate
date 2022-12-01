\name{cybrIDAlleles}
\alias{cybrIDAlleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combine parental alleles with GATK table
%%  ~~function to do ... ~~
}
\description{Function to combine output of cybrConvertParentalAlleles() and cybrQualityFilter() or cyberInputGATKTable().
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrIDAlleles(BSAdfstart, Parentdf, yeast = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSAdfstart}{Dataframe output of either cyberQualityFilter() or cyberInputGATKTable(), required.
%%     ~~Describe \code{x} here~~
}
  \item{Parentdf}{Dataframe output of cyberConvertParentalAlleles().
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{Boolean, if TRUE converts chromosomes to factor, if FALSE does nothing to chromosome names.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Dataframe }{Adds a column to input which contains the SNP ID labeled Allele}
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
