% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanShiftWrapper.R
\name{MeanShiftCluster}
\alias{MeanShiftCluster}
\title{Mean Shift Clustering}
\usage{
MeanShiftCluster(data, ...)
}
\arguments{
\item{\dots}{controls the mean shift}

\item{mydata}{A Data set}
}
\value{
A list of cluster labels and a R object of class "MeanShiftCluster {cluster}"
}
\description{
This function partitions (clustering) of the data into  clusters using the mean shift clustering algorithm.
}
\examples{
library(datasets)
data(iris)

rndSamples <- sample(nrow(iris),100)
trainData <- iris[rndSamples,]
testData <- iris[-rndSamples,]

cls <- MeanShiftCluster(trainData[,1:4])
}
