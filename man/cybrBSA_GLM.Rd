\name{cybrBSA_GLM}
\alias{cybrBSA_GLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to run logistic regression on BSA data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrBSA_GLM(lrP, chr = "II", formula = "PAllele ~ Bulk*Parent",
                           resultscol = c("Intercept", "Bulk", "Parent", "Interaction"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrP}{Input dataframe, containing CHROM, POS, and columns for allele counts and attributes to differentiate by (ie Bulk, Parental Background, or Replicate)
%%     ~~Describe \code{x} here~~
}
  \item{chr}{Character of SINGLE chromosome to include; can be arrayed for all chromosomes individually.
%%     ~~Describe \code{x} here~~
}
  \item{formula}{Input character string for GLM, formatted counts ~ attribute. Attributes and count column must all be column names in input.
%%     ~~Describe \code{x} here~~
}
  \item{resultscol}{Output names as character vector; must include Intercept and each attribute in order. If using an interaction, also include interaction term(s).
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{dataframe }{Output data frame includes CHROM, POS, and columns for each Z-score as defined in input: resultscol vector}
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
Output from prior functions containing Dataset as a column can be used directly as input (in both formula and Results Column) but will only provide the Z score for each dataset at every point. This is useful if running a single BSA of selected/unselected populations.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
