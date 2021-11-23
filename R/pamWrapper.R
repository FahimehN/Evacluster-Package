#' Partitioning Around Medoids (PAM) Clustering
#'
#' This function partitions (clustering) of the data into k clusters "around medoids".
#' In contrast to the k-means algorithm, this clustering methods chooses actual data
#' points as centers
#'
#' @param data A Data set
#' @param k The number of Clusters
#' @return A list of cluster labels and a R object of class "pam {cluster}"
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- pamCluster(trainData[,1:4],3)
#' @export
pamCluster <- function(data=NULL,...)
{
  pm <- cluster::pam(data,...);

  result <- list(classification = pm$clustering,pam = pm);
  class(result) <- "pamCluster"
  return(result);
}

#' pamCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param cls A returned object of pamCluster function
#' @param testData A data set of new samples
#' @return A list of cluster labels
#' @examples
#' pre <- predict(cls,testData[,1:4])
#'
#' @export
predict.pamCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  class <- FRESA.CAD::nearestCentroid(testData,object$pam$medoids)
  result <- list(classification=class)
  return(result);
}
