\name{clac.tumorarray.R}
\alias{clac.tumorarray.R}
\title{A function to do CLAC(Cluster Along Chromosome) Analysis}
\description{A function to do CLAC(Cluster Along Chromosome) Analysis}
   
\usage{
clac.tumorarray.R(NormalResult, tumorarrayIndex)}

\arguments{
    \item{NormalResult}{result object from function \code{clac.preparenormal.R} in package \code{clac}.}
    \item{tumorarrayIndex}{a vector serving as a subscript index to take out a subset of the disease arrays for analysis}
}
        
\details{
    \code{clac.tumorarray.R} do CLAC analysis on selected arrays.
}

\value{
    A list with components
    \item{RegionMean}{data matrix reporting the CLAC result. It has the same dimension as the original input data matrix \code{CANCER} of function \code{clac.preparenormal.R}. For each entry of the matrix, 0 is reported if this genes/clone is not being called to have any copy number alteration, otherwise the amplificaiton/deletion region mean of the genes/clone is reported.}
    \item{fdr}{numeric vector reporting the result FDR for each target disease array.}
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

