\name{cybrConvertParentalAlleles}
\alias{cybrConvertParentalAlleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify parental alleles from VCF Table(s)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrConvertParentalAlleles(
    ParentFiles = c("Wine_VCF.txt", "Oak_VCF.txt"),
    Parents = gsub("_VCF.txt","", ParentFiles),
    Truncate = TRUE, yeast = TRUE
    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ParentFiles}{Character vector of text files to call; must be in VCF format .txt.
%%     ~~Describe \code{x} here~~
}
  \item{Parents}{Character vector of shorthand parental names, default pulls from filenames
%%     ~~Describe \code{x} here~~
}
  \item{Truncate}{Boolean, if TRUE removes positions which do not align to one or the other parental genome; if FALSE keeps all positions. Truncate recommended for downstream functions.
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{Boolean, if TRUE uses the chromosome names for yeast, if FALSE does nothing
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
\item{Dataframe }{Locations of all SNPs corresponding to included parents}
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

