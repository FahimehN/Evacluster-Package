% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getConsensusCluster.R
\name{getConsensusCluster}
\alias{getConsensusCluster}
\title{Consensus Clustering Results}
\usage{
getConsensusCluster(object, who = "training", thr = seq(0.8, 0.3, -0.1))
}
\arguments{
\item{object}{A object of "clusterStability" function result}

\item{who}{This value shows the consensus clustering result of training and testing sets. If who="training" for training set, otherwise other sets.}

\item{thr}{This is the seq function with three arguments that are: initial value, final value, and increment (or decrement for a declining sequence). This produces ascending or descending sequences.}
}
\value{
A list of samples' labels with same connectivity.
}
\description{
This function gets the labels of the subjects that share the same connectivity.
}
\examples{
\donttest{
library("mlbench")
data(Sonar)

Sonar$Class <- as.numeric(Sonar$Class)
Sonar$Class[Sonar$Class == 1] <- 0 
Sonar$Class[Sonar$Class == 2] <- 1

ClustStab <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
                              n_components = 3,featureselection="yes", outcome="Class",
                              fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=3)

clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.30,-0.1))
}
}
