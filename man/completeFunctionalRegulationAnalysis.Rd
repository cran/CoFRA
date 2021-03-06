% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoFRA.R
\name{completeFunctionalRegulationAnalysis}
\alias{completeFunctionalRegulationAnalysis}
\title{This function performs complete functional regulation analysis}
\usage{
completeFunctionalRegulationAnalysis(dfPro, func, Fac, dfComp,
  NbackGround = length(unique(unlist(func))), DataExtract = "FDR",
  minCounts = 10, Test = "t.test", no_cores = -1)
}
\arguments{
\item{dfPro}{data frame with quantitative values}

\item{func}{data frame defining which gene ontology to use "BP","CC","MF"}

\item{Fac}{factor describing the sample groups}

\item{dfComp}{df containing the comparisons to perform}

\item{NbackGround}{integer number of total proteins}

\item{DataExtract}{string which P value correction to use}

\item{minCounts}{integer minimum number of matching genes for functional category}

\item{Test}{"t.test" or "wilcox.test"}

\item{no_cores}{=-1 (no parelle execution) =0 (number of availble cores -1) >0 (use number of cores)}
}
\description{
This function performs complete functional regulation analysis
}
\examples{
library(CoFRA)
data(iBAQ)
Fac=factor(c("MCCTT","MCCTT","MCCTT","MCCT","MCCT","MCCT","MC","MC","MC","MCT","MCT","MCT",
"MTT","MTT","MTT","MT","MT","MT","sN","sN","sN","sNT","sNT","sNT","iN","iN","iN","iNT","iNT","iNT"))
dfComp=data.frame(Con=c("MCCT","MT","MC","iN","sN","AllC,MCCT,MT,MC,iN,sN"),Tre=c("MCCTT","MTT",
"MCT","iNT","sNT","AllT,MCCTT,MTT,MCT,iNT,sNT"))
Func=CoFRA::getFunctionalCategories("CC")
head(str(Func))
CC1=CoFRA::completeFunctionalRegulationAnalysis(iBAQ,Func[100:200],Fac,dfComp,NbackGround=142140) 
}
\keyword{heatmap}

