#' K-means Clustering
#'
#' This function classifies unlabeled data by grouping them by features,
#' rather than pre-defined categories. It splits the data into K different
#' clusters and describes the location of the center of each cluster. Then,
#' a new data point can be assigned a cluster (class) based on the closed
#' center of mass.
#'
#' @param data A Data set
#' @param \dots center: The number of centers
#' @return A list of cluster labels and a R object of class "kmeans"
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- kmeansCluster(trainData[,1:4],3)
#' @export
kmeansCluster <- function(data=NULL,...)
{
  km <- stats::kmeans(data,...);

  result <- list(classification = km$cluster,kmeans = km);
  class(result) <- "kmeansCluster"
  return(result);
}

#' kmeansCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of kmeansCluster function
#' @param \dots New samples set
#' @return A list of cluster labels
#'
#' @export
predict.kmeansCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  class <- FRESA.CAD::nearestCentroid(testData,object$kmeans$centers)
  result <- list(classification=class)
  return(result);
}
