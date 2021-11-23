#' @export
gmmCluster <- function(data=NULL,gaussian_comps=NULL,dist_mode=NULL,seed_mode=NULL,km_iter,em_iter,verbose,...)
{

  gmm <- ClusterR::GMM(data, gaussian_comps, dist_mode, seed_mode, km_iter,em_iter, verbose,...)

  result <- list(centroids = gmm$centroids,covariance = gmm$covariance_matrices, weights=gmm$weights,gmm = gmm);
   class(result) <- "gmmCluster"
  return(result);
}

predict.gmmCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  class <- ClusterR::predict_GMM(testData, object$centroids, object$covariance, object$weights)
  result <- list(classification=class$cluster_labels+1)
  return(result);
}
