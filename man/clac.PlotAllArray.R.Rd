\name{clac.PlotAllArray.R}
\alias{clac.PlotAllArray.R}
\title{ A function to plot multiple arrays for CLAC result.}
\description{A function to plot multiple arrays for CLAC result.}
   
\usage{
clac.PlotAllArray.R(NormalResult, clac.result, ArrayName=NULL,centromere=NULL)
}
\arguments{
    \item{NormalResult}{result object of function \code{clac.preparenormal.R}}
    \item{clac.result}{result object of function clac.tumorarray.R}
    \item{ArrayName}{string vector, providing the name of each disease array in \code{clac.result}. The length should be the same as the total number of disease arrays.}
    \item{centromere}{numeric vector specifying the centromere positions. If missing, the default centromere value of human genome will be used.}
 }  
    
\details{
    \code{clac.PlotAllArray.R} plot multiple array for CLAC result.}

\value{
    No return value.
   }

\references{
P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani, ¡°A method for calling gains and losses in array CGH data¡±, Biostatistics (accepted for publication 4/5/2004), available at http://www-stat.stanford.edu/~wp57/CGH-Miner/
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
