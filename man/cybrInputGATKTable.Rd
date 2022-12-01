\name{cybrInputGATKTable}
\alias{cybrInputGATKTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to convert VCF datatable into R dataframe
%%  ~~function to do ... ~~
}
\description{
Converts the gatk output into a data frame, renaming the AD, PL, GQ, and DP fields for each sample in the file. Will create a column for dataset based on each unique filename, and will split the AD and PL by alternate and reference alleles. Genome quality can be filtered by using function cybrQualityFilter() on the output. See Details for GATK script to create input.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrInputGATKTable(rawData, yeast = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawData}{Output of gatk VariantsToTable, required. See Details for how to create table. Columns required: CHROM (chromosome), POS (Position), REF (reference allele frequency), ALT (alternate allele frequency), and PER SAMPLE AD (overall allele depth), DP (?), GQ (genome quality), PL (?).
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{If TRUE, converts chromosome names to Roman Numerals; if FALSE, does no conversion
%%     ~~Describe \code{x} here~~
}
}
\details{

gatk VariantsToTable \
     -V ${myfile} \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O ${myfile}.output.table

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Dataframe}{Output contains columns CHROM, POS, Dataset (unique names), Allele (ALT/REF), AD (overall allele depth), DP (?), GQ (genome quality), and PL (Phred-scaled likelihoods of the genotypes)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res, 20:1297-303. DOI: 10.1101/gr.107524.110.
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
