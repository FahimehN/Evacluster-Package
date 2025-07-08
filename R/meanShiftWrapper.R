#' Mean Shift Clustering
#'
#' This function partitions (clustering) of the data into  clusters using the mean shift clustering algorithm.
#' 
#' @param mydata A Data set
#' @param \dots controls the mean shift
#' @return A list of cluster labels and a R object of class "MeanShiftCluster {cluster}"
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- MeanShiftCluster(trainData[,1:4])
#' @export
MeanShiftCluster <- function(data,...)
{
  if (!requireNamespace("meanShiftR", quietly = TRUE)) {
    install.packages("meanShiftR", dependencies = TRUE)
  } 
  if (!requireNamespace("FRESA.CAD", quietly = TRUE)) {
    install.packages("FRESA.CAD", dependencies = TRUE)
  } 
  if (!requireNamespace("MASS", quietly = TRUE)) {
    install.packages("MASS", dependencies = TRUE)
  } 

  data <- as.matrix(data);
  parameters <- list(...)
#  print(c(parameters,length(parameters)))
  if (length(parameters)==0)
  {
#    cat("W")
    cluster <- meanShiftR::meanShift(data,
                         nNeighbors=round(0.40*nrow(data)),
                         iterations=20,
                         alpha=0.0,
                         epsilon=1.0e-8,
                         epsilonCluster=1.0e-2,
                         bandwidth=rep(0.35,NCOL(data)));
  }
  else
  {
    cluster <- meanShiftR::meanShift(data,
                                     ...);
    
  }
  
  numlabesl <- unique(cluster$assignment)
  meanV = list()
  covM = list()
  lbt <- 0;
  smallestCluster <- max(c(ncol(data) + (ncol(data)^2)/2,0.01*nrow(data)));
  for (lb in numlabesl)
  {
    dtlab = data[cluster$assignment==lb,];
#    cat(length(dtlab),",",nrow(dtlab),"\n")
    if (length(dtlab) > ncol(data))
    {
      if (nrow(dtlab) > smallestCluster)
      {
        lbt <- lbt+1;
#        meanV[[lbt]] = apply(dtlab,2,mean);
#        covM[[lbt]] = cov(dtlab);
        mve_fit <- MASS::cov.rob(dtlab,method = "classical")
        meanV[[lbt]] = mve_fit$center;
        covM[[lbt]] = mve_fit$cov;
      }
    }
  }
#  cat(smallestCluster,":number of clsters:(",numlabesl,",",lbt,")\n");
  result <- list(cluster = cluster,meanV=meanV,covM=covM)
  class(result) <- "MeanShiftCluster"
  return(result)
}


#' MeanShiftCluster prediction function
#'
#' This function predicts the labels of the cluster for new data based on
#' cluster labels of the training set.
#'
#' @param object A returned object of MeanShiftCluster function
#' @param \dots New samples set
#' @return A list of cluster labels
#'
#' @export
predict.MeanShiftCluster <- function(clusResult,...)
{
  parameters <- list(...)
  data <- as.matrix(parameters[[1]]);
  testlabesl <- FRESA.CAD::nearestCentroid(data,clusResult$meanV,clusResult$covM,p.threshold = 0);
  names(testlabesl) <- rownames(data);
  result <- list(classification=testlabesl)
  return(result)
}
