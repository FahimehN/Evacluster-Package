% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterStability.r
\name{clusterStability}
\alias{clusterStability}
\title{clustering stability function}
\usage{
clusterStability(
  data = NULL,
  clustermethod = MeanShiftCluster,
  dimenreducmethod = c("Auto", "none", "PCA", "tsne", "UMAP"),
  n_components = 3,
  perplexity = 25,
  max_iter = 1000,
  k_neighbor = 3,
  featureselection = NULL,
  outcome = NULL,
  fs.pvalue = 0.05,
  randomTests = 20,
  trainFraction = 0.5,
  pac.thr = 0.1,
  plotClustering = FALSE,
  ...
)
}
\arguments{
\item{data}{A Data set}

\item{clustermethod}{The clustering method. This can be one of "Mclust","pamCluster","kmeansCluster", "hierarchicalCluster",and "FuzzyCluster".}

\item{dimenreducmethod}{The dimensionality reduction method. This must be one of "Auto", "none","UMAP","tSNE", and "PCA".}

\item{n_components}{The dimension of the space that data embed into. It can be set to any integer value in the range of 2 to 100.}

\item{perplexity}{The Perplexity parameter that determines the optimal number of neighbors in tSNE method.(it is only used in the tSNE reduction method)}

\item{max_iter}{The maximum number of iterations for performing tSNE reduction method.}

\item{k_neighbor}{The k_neighbor is used for computing the means of #neighbors with min distance (#Neighbor=sqrt(#Samples/k) for performing an embedding of new data using an existing embedding in the tSNE method.}

\item{featureselection}{This parameter determines whether feature selection is applied before clustering data or not. if used, it should be "yes", otherwisw "no".}

\item{outcome}{The outcome feature is used for feature selection.}

\item{fs.pvalue}{The threshold pvalue used for feature selection process. The default value is 0.05.}

\item{randomTests}{The number of iterations of the clustering process for computing the cluster stability.}

\item{trainFraction}{This parameter determines the ratio of training data. The default value is 0.5.}

\item{pac.thr}{The pac.thr is the threshold to use for computing the proportion of ambiguous clustering (PAC) score. It is as the fraction of sample pairs with consensus indices falling in the interval.The default value is 0.1.}

\item{plotClustering}{if TRUE the class-labeled scatter plot of the first two dimensions will be shown.}

\item{\dots}{Additional arguments passed to the clustering algorithm.}
}
\value{
A list with the following elements:
\itemize{
  \item randIndex - A vector of the Rand Index that computes a similarity measure between two clusterings. 
  \item jaccIndex - A vector of jaccard Index that measures how frequently pairs of items are joined together in two clustering data sets.
  \item randomSamples - A vector with indexes of selected samples for training in each iteration.
  \item clusterLabels - A vector with clusters' labels in all iterations. jaccardpoint
  \item jaccardpoint - The corresponding Jaccard index for each data point of testing set
  \item averageNumberofClusters - The mean Number of Clusters.
  \item testConsesus - A vector of consensus clustering results of testing set.
  \item trainRandIndex - A vector of the Rand Index for training set.
  \item trainJaccIndex - A vector of the jaccard Index for training set.
  \item trainJaccardpoint - The corresponding Jaccard index for each data point of training set.
  \item PAC - The proportion of ambiguous clustering (PAC) score.
  \item dataConcensus - A vector of consensus clustering results of training set.
}
}
\description{
This function computes the stability of clustering that helps to select the best number of clusters.
Feature selection and dimensionality reduction methods can be used before clustering 
the data.
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


ClustStab <- clusterStability(data=Sonar, clustermethod=pamCluster, dimenreducmethod="tSNE",
                              n_components = 3, perplexity=10,max_iter=100,k_neighbor=2,
                              featureselection="yes", outcome="Class",fs.pvalue = 0.05,
                              randomTests = 100,trainFraction = 0.7,k=3)


ClustStab <- clusterStability(data=Sonar, clustermethod=hierarchicalCluster, 
                              dimenreducmethod="PCA", n_components = 3,featureselection="no",
                              randomTests = 100,trainFraction = 0.7,distmethod="euclidean",
                              clusters=3)

}
}
