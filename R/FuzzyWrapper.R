#' Fuzzy C-means Clustering Algorithm
#'
#' This function works by assigning membership to each data point corresponding
#' to each cluster center based on the distance between the cluster center and
#' the data point. A data object is the member of all clusters with varying
#' degrees of fuzzy membership between 0 and 1.
#'
#' @param data A Data set
#' @param \dots k: The number of Clusters
#' @return A list of cluster labels and a R object of class "fcm {ppclust}"
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- FuzzyCluster(trainData[,1:4],3)
#' @export
FuzzyCluster <- function(data=NULL,...)
{
  parameters <- list(...);
  k <- parameters[[1]];

  v <- inaparc::kmpp(data,k)$v
  u <- inaparc::imembrand(nrow(data),k)$u
  fc <- ppclust::fcm(data,centers=v, memberships=u,...)
  
  result <- list(classification = fc$cluster,fuzzy = fc);
  class(result) <- "FuzzyCluster"
  return(result);
}

#' FuzzyCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of FuzzyCluster function
#' @param \dots New samples set
#' @return A list of cluster labels
#'
#' @export
predict.FuzzyCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  class <- FRESA.CAD::nearestCentroid(testData,object$fuzzy$v0)
  result <- list(classification=class)
  return(result);
}
