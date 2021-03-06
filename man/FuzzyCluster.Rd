% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuzzyWrapper.R
\name{FuzzyCluster}
\alias{FuzzyCluster}
\title{Fuzzy C-means Clustering Algorithm}
\usage{
FuzzyCluster(data = NULL, ...)
}
\arguments{
\item{data}{A Data set}

\item{\dots}{k: The number of Clusters}
}
\value{
A list of cluster labels and a R object of class "fcm {ppclust}"
}
\description{
This function works by assigning membership to each data point corresponding
to each cluster center based on the distance between the cluster center and
the data point. A data object is the member of all clusters with varying
degrees of fuzzy membership between 0 and 1.
}
\examples{
library(datasets)
data(iris)

rndSamples <- sample(nrow(iris),100)
trainData <- iris[rndSamples,]
testData <- iris[-rndSamples,]

cls <- FuzzyCluster(trainData[,1:4],3)
}
