\name{clac.PlotConsensus.R}
\alias{clac.PlotConsensus.R}
\title{ A function to make consensus plot for CLAC result.}
\description{A function to make consensus plot for CLAC result.}
   
\usage{
clac.PlotConsensus.R(clac.result, chromosome.number, nucleotide.position, sample.index, centromere=NULL,graylevel=0.9)
}

\arguments{
    \item{clac.result}{Result object of function clac.tumorarray.R}
    \item{chromosome.number}{numeric vector . Length should be the same as the row number of \code{NORMAL}. It's the chromosome number of each gene/clone.}
    \item{nucleotide.position}{numeric vector. Length should be the same as the row number of \code{NORMAL}. It's the nucleotide position of each gene/clone.}  
    \item{sample.index}{integer vector. It's the index of the arrays to be included in the consensus plot.}
    \item{centromere}{numeric vector specifying the centromere positions. If missing, the default centromere value of human genome will be used.}
    \item{graylevel}{numeric value between 0 and 1, specifying the gray level for the background scale lines. 0 for black and 1 for white.}
 }  
    
\details{
    \code{clac.PlotConsensus.R} make consensus plot for CLAC result.
 }

\value{
    No return value.
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

