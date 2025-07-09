#' Consensus Clustering Results
#'
#' This function gets the labels of the subjects that share the same connectivity.
#'
#' @param object A object of "clusterStability" function result
#' @param who This value shows the consensus clustering result of training and testing sets. If who="training" for training set, otherwise other sets. 
#' @param thr This is the seq function with three arguments that are: initial value, final value, and increment (or decrement for a declining sequence). This produces ascending or descending sequences.
#' @return A list of samples' labels with same connectivity. 
#' Additional attributes include:
#' \describe{
#'   \item{Mean_in}{The average Co-association in the cluster}
#'   \item{Mean_out}{The average co-association outside the cluster}
#'   \item{SD_in}{The standard deviation of the  Co-association in the cluster}
#'   \item{SD_out}{The standard deviation of the  Co-association outside the cluster}
#'   \item{Quality}{The quality of each cluster}
#' }
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
  #  cat(thr,"Here\n")
  for (lthr in thr)
  {
    totlabeled <- sum(classID > 0);
    #    print(totlabeled)
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
            #            print(c(SconA,SconB))
            
            if ( (SconB > 0.01*npoints) || (SconA > 0.01*npoints) )
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
      cat(minLables,":",sprintf("%5.3f",lthr),": ",totlabeled,": \n")
    }
  }
  classID <- classID[orgnames];
  theClasses <- table(classID)
  avgindx <- numeric(length(theClasses))
  stdidx <- numeric(length(theClasses))
  navgindx <- numeric(length(theClasses))
  nstdidx <- numeric(length(theClasses))
  ix <- 0;
  for (cid in names(theClasses))
  {
    ix <- ix + 1
    whosub <- names(classID)[classID==cid];
    notsub <- names(classID)[classID!=cid];
    avgindx[ix] <- mean(concensusMat[whosub,whosub]);
    stdidx[ix] <- sd(concensusMat[whosub,whosub]);
    navgindx[ix] <- mean(concensusMat[whosub,notsub]);
    nstdidx[ix] <- sd(concensusMat[whosub,notsub]);
  }
  quality <- rep(1.0,length(theClasses))-sqrt((stdidx^2+nstdidx^2))/(avgindx-navgindx);
  names(avgindx) <- names(theClasses)
  names(stdidx) <- names(theClasses)
  names(navgindx) <- names(theClasses)
  names(nstdidx) <- names(theClasses)
  names(quality) <- names(theClasses)
  attr(classID,"Mean_in") <- avgindx
  attr(classID,"SD_in") <- stdidx
  attr(classID,"Mean_out") <- navgindx
  attr(classID,"SD_out") <- nstdidx
  attr(classID,"Quality") <- quality
  return (classID);
}
