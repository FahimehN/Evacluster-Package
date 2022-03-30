#' Consensus Clustering Results
#'
#' This function gets the labels of the subjects that share the same connectivity.
#'
#' @param object A object of "clusterStability" function result
#' @param who This value shows the consensus clustering result of training and testing sets. If who="training" for training set, otherwise other sets. 
#' @param thr This is the seq function with three arguments that are: initial value, final value, and increment (or decrement for a declining sequence). This produces ascending or descending sequences.
#' @return A list of samples' labels with same connectivity.
#' @examples
#' \donttest{
#' library("mlbench")
#' data(Sonar)
#' 
#' Sonar$Class <- as.numeric(Sonar$Class)
#' Sonar$Class[Sonar$Class == 1] <- 0 
#' Sonar$Class[Sonar$Class == 2] <- 1
#'
#' ClustStab <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
#'                               n_components = 3,featureselection="yes", outcome="Class",
#'                               fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=3)
#'
#' clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.30,-0.1))
#' }
#' @export
getConsensusCluster <- function(object,who="training",thr=seq(0.80,0.30,-0.1))
{
  
  orgnames <-  rownames(object$dataConcensus);
  if (who != "training")
  {
    orgnames <-  rownames(object$testConsesus);
    pointJaccard <- object$jaccardpoint;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$testConsesus[order(-pointJaccard),]
  }
  else
  {
    pointJaccard <- (object$trainJaccardpoint+object$jaccardpoint)/2.0;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$dataConcensus[order(-pointJaccard),]
  }
  concensusMat <- concensusMat[,order(-pointJaccard)]
  classID <- numeric(nrow(concensusMat));
  names(classID) <-  rownames(concensusMat);
  pointJaccard <- pointJaccard[order(-pointJaccard)];
  npoints <- length(pointJaccard)
  label <- 1;
  for (lthr in thr)
  {
    totlabeled <- sum(classID > 0);
    if (totlabeled < npoints)
    {
      added <- 1;
      while (added > 0)
      {
        added <- 0;
        for (i in 1:npoints)
        {
          if (classID[i] == 0)
          {
            minLables <- label;
            wcon <- concensusMat[i,];
            consensA <- (wcon > lthr) & (classID > 0)
            consensB <- (wcon > lthr) & (classID == 0)
            SconA <- sum(pointJaccard[consensA]);
            SconB <- sum(pointJaccard[consensB]) - pointJaccard[i];
            
            if ( (SconB > 0.05*npoints) || (SconA > 0.05*npoints) )
            {
              if (SconB > SconA)
              {
                classID[consensB] <- label;
                added <- 1;
                label <- label + 1;
              }
              else
              {
                if (SconB >= 0)
                {
                  if (SconA > 0)
                  {
                    tb <- table(classID[consensA])
                    minLables <- as.numeric(names(which.max(tb))[1])
                    if (sum(pointJaccard[classID == minLables]) < SconB )
                    {
                      minLables <- label;
                    }
                  }
                  classID[consensB] <- minLables;
                  added <- 1;
                  if (minLables == label)
                  {
                    label <- label + 1;
                  }
                }
              }
            }
          }
        }
      }
      totlabeled <- sum(classID > 0);
      warning(minLables,":",sprintf("%5.3f",lthr),": ",totlabeled,": ")
    }
  }
  classID <- classID[orgnames];
  return (classID);
}
