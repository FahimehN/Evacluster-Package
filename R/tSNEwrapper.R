#' t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' This method is an unsupervised, non-linear technique used for
#' data exploration and visualizing high-dimensional data.This function
#' constructs a low-dimensional embedding of high-dimensional
#' data, distances, or similarities.
#'
#' @param data Data matrix (each row is an observation, each column is a variable)
#' @param dim  Integer number; Output dimensional (default=2)
#' @param perplexity  numeric; Perplexity parameter (should not be bigger than
#'                     3 * perplexity < nrow(X) - 1, default=30)
#' @param max_iter  Integer; Number of iterations (default: 500)
#' @return tsneY: A Matrix containing the new representations for the observation
#'          with selected dimensions by user
#' @examples
#' library("mlbench")
#' data(Sonar)
#' 
#' rndSamples <- sample(nrow(Sonar),150)
#' trainData <- Sonar[rndSamples,]
#' testData <- Sonar[-rndSamples,]
#'
#' tsne_trainData <- tsneReductor(trainData[,1:60],dim = 3,perplexity = 10,max_iter = 1000)
#'
#' @export
tsneReductor <- function(data=NULL,dim=2,perplexity=30,max_iter=500)
{
  Tsne <- Rtsne::Rtsne(unique(data), dims = dim, perplexity=perplexity, max_iter = max_iter);#, verbose=TRUE

  result <- list(Data = unique(data) ,tsneY = Tsne$Y);
  class(result) <- "tsneReductor"
  return(result);
}


#' tsneReductor prediction function
#'
#' This function performs an embedding of new data using an existing embedding.
#'
#' @param object   A returned object of tsneReductor function
#' @param k   The number is used for computing the means of #neighbors with min distance
#'          (#Neighbor=sqrt(#Samples/k).
#' @param \dots New samples set
#' 
#' @return tsneY:An embedding of new data
#' @examples
#' library("mlbench")
#' data(Sonar)
#' 
#' rndSamples <- sample(nrow(Sonar),150)
#' trainData <- Sonar[rndSamples,]
#' testData <- Sonar[-rndSamples,]
#'
#' tsne_trainData <- tsneReductor(trainData[,1:60],dim = 3,perplexity = 10,max_iter = 1000)
#' 
#' tsne_testData <- predict(tsne_trainData,k=3,testData[,1:60])
#' @export

predict.tsneReductor <- function(object,k=NULL,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];

  dist <- proxy::dist(object$Data,testData, method = "euclidean")

  mindist <- apply(dist,2,min)

  minDist <- apply(object$tsneY,2,min)
  maxDist <- apply(object$tsneY,2,max)

  thr_val <- ((maxDist[1] - minDist[1])/nrow(testData))*2

  prdtSNE <- numeric();
  Lthr <- 0;
  Mthr <- 0;
  N <- 0;

  for(i in 1:length(mindist))
  {
    sidx <- which(dist[,i] == mindist[i])

    if(dist[sidx,i] < thr_val) {
      prdtSNE <- rbind(prdtSNE,object$tsneY[sidx,])
      Lthr <- Lthr + 1
    } else {
      N <- sqrt(nrow(object$Data)/k)
      if(N<1){
        N <- 1}
      ndx <- order(dist[,i])[1:N]
      prdtSNE <- rbind(prdtSNE,colMeans(object$tsneY[ndx,]))
      Mthr <- Mthr + 1
    }
  }

  result <- list(tsneY=prdtSNE,Lthr=Lthr,Mthr=Mthr,minDist=minDist,
                         maxDist=maxDist,thr=thr_val,N=N)
  return(result);

}

