\name{CLAC analysis for CGH-Miner}
\alias{clac.from.Excel}
\alias{clac.preparenormal.Excel}
\alias{clac.onetumorarray.Excel}
\alias{clac.finalsummary.Excel}
\alias{clac.PlotOneArray}
\alias{clac.PlotConsensus}

\title{R Functions for Excel Addin Packgage CGH-Miner}

\description{These functions are desinged to be called by Excel Addin package CGH-Miner to perform CLAC analysis.}
   
\usage{
clac.from.Excel(column.name, SampleID, Data, windowsize, FDR, log2trans, centerColumn, chromosomeOption)
clac.preparenormal.Excel(inputdata)
clac.onetumorarray.Excel(NormalResult, i)
clac.finalsummary.Excel(NormalResult, CLAC.result)
clac.PlotOneArray(inputdata, clac.result, Sample) 
clac.PlotConsensus(inputdata, clac.result, SampleIDPick)
}

\arguments{
    \item{column.name}{string vector, provides the column names in the spread sheet.}
    \item{SampleID}{string vector, provides the sample ID for each CGH array.}
    \item{Data}{numeric matrix, contains the data read in from the spread sheet.}
    \item{windowsiz}{enumeric value, specifying the window size to carry out the average smooth.}
    \item{FDR}{numeric value between 0 and 1, specifying the desired fianl FDR for CLAC analysis.}
    \item{log2trans}{boolean variable, if TRUE, a log2 transform will be performed on \code{Data}.}
    \item{centerColumn}{boolean variable, if TRUE, the columns of \code{Data} will be centered to zero.}
    \item{chromosomeOption}{boolean variable. If False, the chromosome arms will be considered seperately. If true, two chormosome arms of one chromosome would be dealed together.}
    \item{inputdata}{retruned object from function \code{clac.from.Excel}.}
    \item{NormalResult}{return object from function \code{clac.preparenormal.Excel}.}
    \item{CLAC.result}{numeric matrix. First row records the estimated FDR for each array; second row records the \code{alpha} value for each array; all left rows records the result of CLAC analysis for each array.}
    \item{Sample}{string variable, is the sample ID of one array.}
    \item{SampleIDPick}{string vector, is a subset of \code{SampleID}.}
}
    
\details{
    These functions are not to be used in RGui interface directly. For \code{clac} functions in R, please see \code{clac.preparenormal.R},
\code{clac.tumorarray.R, clac.PlotSingleArray.R, clac.PlotConsensus.R, clac.PlotAllArray.R}.
}

\references{
P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani, ?A method for calling gains and losses in array CGH data?, Biostatistics (accepted for publication 4/5/2004), available at http://www-stat.stanford.edu/~wp57/CGH-Miner/
}

\keyword{internal}


\author{Pei Wang}



