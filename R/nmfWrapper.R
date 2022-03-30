#' Non-negative matrix factorization (NMF)
#'
#' This function factorizes samples matrix into (usually) two matrices W the cluster centroids
#'  and H the cluster membership,
#' @param data A Data set
#' @param rank Specification of the factorization rank
#'
#' @return A list of cluster labels, a R object of class "nmf" and the centers of the clusters
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- nmfCluster(trainData[,1:4],rank=3)
#' @export
nmfCluster <- function(data=NULL,rank=NULL)
{
  
  nmf <- NMF::nmf(t(data),rank);
  
  W <- NMF::basis(nmf);#gives the cluster centroids
  H <- NMF::coef(nmf); #gives the cluster membership

  clusters <- NMF::predict(nmf)

  result <- list(Data = data,classification = clusters,centers = W,nmf = nmf,H = H,Rank = rank);
  class(result) <- "nmfCluster"
  return(result);
}

#' nmfCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of nmfCluster
#' @param \dots New samples set
#' @return A list of cluster labels
#'
#' @export
predict.nmfCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];

  nmf <- NMF::nmf(t(testData),object$Rank);

  W <- NMF::basis(nmf);
  iW <- MASS::ginv(W)

  H <- iW %*% t(testData)

  class <- apply(H, 2, which.max)
  result <- list(classification=class)
  return(result);
}
