#' hierarchical clustering
#'
#' This function seeks to build a hierarchy of clusters
#'
#' @param data A Data set
#' @param distmethod The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param clusters The number of Clusters
#' @param \dots Additional parameters passed to  hclust function
#' @return A list of cluster labels
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- hierarchicalCluster(trainData[,1:4],distmethod="euclidean",clusters=3)
#' @export
hierarchicalCluster <- function(data=NULL,distmethod=NULL,clusters=NULL,...)
{
  distance= stats::dist(data,distmethod);
  hc0= stats::hclust(d = distance,...);
  k <- clusters;
  hc <- stats::cutree(hc0,k);
  
  result <- list(Data = data ,classification = hc);
  class(result) <- "hierarchicalCluster"
  return(result);
}

#' hierarchicalCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of hierarchicalCluster function
#' @param \dots New samples set
#' @return A list of cluster labels
#'
#' @export
predict.hierarchicalCluster <- function(object,...)
{
  parameters <- list(...);
  if(is.null(parameters$kn))
  {
    kn <- as.integer(sqrt(nrow(object$Data))+0.5)
  }
  else{
    kn <- parameters$kn
  }
  testData <- parameters[[1]];
  class <- class::knn(object$Data,testData ,factor(object$classification),kn)
  result <- list(classification=as.integer(as.character(class)))
  return(result);
}
