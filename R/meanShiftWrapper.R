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
  

  data <- as.matrix(data);
  parameters <- list(...)
  cluster <- meanShift(data,
                       nNeighbors=round(0.40*nrow(data)),
                       iterations=20,
                       alpha=0.0,
                       epsilon=1.0e-8,
                       epsilonCluster=1.0e-2,
                       bandwidth=rep(0.35,NCOL(data)));
  
  numlabesl <- unique(cluster$assignment)
  meanV = list()
  covM = list()
  lbt <- 0;
  smallestCluster <- max(c(ncol(mydata) + (ncol(mydata)^2)/2,0.01*nrow(mydata)));
  for (lb in numlabesl)
  {
    dtlab = mydata[cluster$assignment==lb,];
    if (length(dtlab) > 3)
    {
      if (nrow(dtlab) > smallestCluster)
      {
        lbt <- lbt+1;
        meanV[[lbt]] = apply(dtlab,2,mean);
        covM[[lbt]] = cov(dtlab);
      }
    }
  }
  cat(smallestCluster,":number of clsters:(",numlabesl,",",lbt,")\n");
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
  testlabesl <- nearestCentroid(data,clusResult$meanV,clusResult$covM,p.threshold = 0);
  names(testlabesl) <- rownames(data);
  result <- list(classification=testlabesl)
  return(result)
}
