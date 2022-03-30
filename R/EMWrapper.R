#' Expectation Maximization Clustering
#'
#' This function perform EM algorithm for model-based clustering of finite
#' mixture multivariate Gaussian distribution.The general purpose of clustering
#' is to detect clusters of data and to assign the data to the clusters.
#'
#' @param data A Data set
#' @param \dots k: The number of Clusters
#' @return A list of cluster labels and a returned object from init.EM
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' clsut <- EMCluster(trainData[,1:4],3)
#' @export
EMCluster <- function(data=NULL,...)
{
  em <- EMCluster::exhaust.EM(data,...)
  
  result <- list(classification = as.integer(em$class),EM = em);
  class(result) <- "EMCluster"
  return(result);
}

#' EMCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of EMCluster
#' @param \dots New sample set
#' @return A list of cluster labels
#' 
#' @export
predict.EMCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  pr <- EMCluster::assign.class(testData,object$EM,return.all = FALSE)
  result <- list(classification=pr$class)
  return(result);
}
