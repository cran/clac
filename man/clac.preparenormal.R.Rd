\name{clac.preparenormal.R}
\alias{clac.preparenormal.R}
\title{ A function to prepare normal reference in CLAC(Cluster Along Chromosome) Analysis}
\description{
CLAC is a method for calling gains and losses in CGH array data.
 This function is to prepare normal references in CLAC(Cluster Along Chromosome) Analysis}
   
\usage{
clac.preparenormal.R(CANCER, NORMAL, Normal.Type, chromosome.number, nucleotide.position, windowsize=5, targetFDR=0.01, chromosomeOption=FALSE, centromere=NULL)
}

\arguments{
    \item{CANCER}{data matrix. It's the result of a group of CGH experiments, which are the target disease arrays to analysis. Each column corresponds to one sample (one target array), and each row corresponds to one gene/clone. The (i, j) entry should be the log fluorescence ratio of the ith gene/clone in the jth sample. Missing value should be coded as either NA or 999.}
    \item{NORMAL}{data matrix. It's also the result of a group of CGH experiment. But these results are from normal reference arrays. Again, each column corresponds to one sample, and each row corresponds to one gene/clone. Missing value should be coded as either NA or 999.}
    \item{Normal.Type}{a vector specifying the normal array type. Length should be the same as the column number of \code{NORMAL}. Code 0 for normal reference arrays from the same gender hybridization, while 1 for arrays from opposite gender hybridization. }
    \item{chromosome.number}{numeric vector . Length should be the same as the row number of \code{NORMAL}. It's the chromosome number of each gene/clone.}
    \item{nucleotide.position}{numeric vector. Length should be the same as the row number of \code{NORMAL}. It's the nucleotide position of each gene/clone.}  
    \item{windowsize}{numeric value, specifying the window size to carry out the average smooth.}
    \item{targetFDR}{numeric value between 0 and 1, specifying the desired fianl FDR for CLAC analysis.}
    \item{chromosomeOption}{a boolean variable. If False, the chromosome arms will be considered seperately. If true, two chormosome arms of one chromosome would be dealed together.}
    \item{centromere}{numeric vector specifying the centromere positions. If missing, the default centromere value of human genome will be used.}
}  
    
\details{
    \code{clac.preparenormal.R} builds cluster trees on normal reference arrays, the result is reported to the next step of CLAC analysis.}

\value{
    A list with components
    \item{normal.result}{An object containing the information about the normal reference arrays.}
    \item{CANCER.sm}{The result of average smooth for data matrix CANCER}
    \item{NORMAL.sm}{The result of average smooth for data matrix NORMAL}
   }

\references{
P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani, ?A method for calling gains and losses in array CGH data?, Biostatistics (accepted for publication 4/5/2004), available at http://www-stat.stanford.edu/~wp57/CGH-Miner/
}

\author{ Pei Wang}

\keyword{cluster}

\examples{
library(clac)
data(BACarray)
attach(BACarray)

############ prepare the normal reference arrays
NormalResult<-clac.preparenormal.R(DiseaseArray, NormalArray, Normal.Type=rep(0,3), chromosome.number=chromosome, nucleotide.position=nucposition, windowsize=5, targetFDR=0.01, chromosomeOption=FALSE)

############ clac on selected tumor arrays
clac.result<-clac.tumorarray.R(NormalResult, tumorarrayIndex=1:4)

############ Plot for the first arrays
i<-1
clac.PlotSingleArray.R(i, NormalResult, clac.result)
title(main=paste("CLAC Plot for the ", i ,"th BAC array; FDR=", round(clac.result$fdr[i],3), sep=""))

############ consensus plot
clac.PlotConsensus.R(clac.result, chromosome, nucposition,  1:4)
title(main="Consensus Plot for 4 BAC arrays")

############ Plot all arrays
clac.PlotAllArray.R(NormalResult, clac.result)
title(main="Plot for all 4 arrays")
}
