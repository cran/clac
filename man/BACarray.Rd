\name{BACarray}
\alias{BACarray}
\title{Example CGH BAC Array Data}
\description{A list containing an example CGH BAC array data for package \code{clac}}

\format{
BACarray is a list of four components: 
\tabular{rl}{
NormalArray \tab  a numeric matrix consisting of 2270 rows and 3 columns \cr
DiseaseArray \tab a numeric matrix consisting of 2270 rows and 4 columns \cr 
chromosome  \tab  a numeric vector of length 2270\cr
nucposition  \tab a numeric vector of length 2270
}
}

\details{
In \code{NormalArray} and \code{DiseaseArray}, each column corresponds to
one CGH array, and each row corresponds to one gene/clone. The value of each entry is the log fluorescence
ratio resulted from the CGH experiment. The order of the genes/clones in the rows is the same as the order of the
genes/clones on the genome. \code{chromosome} and \code{nucposition} provide chromosome number and nucleotide 
position for each gene/clone.}

\references{
P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani, ¡°A method for calling gains and losses in array CGH data¡±, Biostatistics (accepted for publication 4/5/2004), available at http://www-stat.stanford.edu/~wp57/CGH-Miner/
}

\usage{data(BACarray)}

\keyword{datasets}

