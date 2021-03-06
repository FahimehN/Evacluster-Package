% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmfWrapper.R
\name{predict.nmfCluster}
\alias{predict.nmfCluster}
\title{nmfCluster prediction function}
\usage{
\method{predict}{nmfCluster}(object, ...)
}
\arguments{
\item{object}{A returned object of nmfCluster}

\item{\dots}{New samples set}
}
\value{
A list of cluster labels
}
\description{
This function predicts the labels of the cluster for new data based on
cluster labels of the training set.
}
