#' mergeConsensusClusters
#'
#' This function merge clusters that share consensus.
#'
#' @param object A object of "getConsensusCluster" function result
#' @param Zthr below this z-value two clusters will be merged
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
#' clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.0,-0.1))
#' NewclusterLabels <- getConsensusCluster(clusterLabels)
#' }
#' @export
mergeConsensusClusters <- function(object,Zthr=2.0)
{
  
  concensusMat <- attr(object,"concensusMat");
  pointJaccard <- attr(object,"pointJaccard");
  orgnames <-  rownames(concensusMat);

  classID <- object[orgnames];
  theClasses <- table(classID)
  oclases <- names(theClasses)
  sclases <- oclases; 
  for (cid in names(theClasses))
  {
    if (sum(cid %in% sclases) > 0)
    {
      oclases <- oclases[!(oclases %in% cid)];
      for (oid in oclases)
      {
        if (sum(oid %in% sclases) > 0)
        {
          whosub <- names(classID)[classID==cid];
          otrsub <- names(classID)[classID==oid];
          meanwho <- mean(concensusMat[whosub,whosub]);
          stdwho <- sd(concensusMat[whosub,whosub]);
          meanotr <- mean(concensusMat[otrsub,otrsub]);
          stdotr <- sd(concensusMat[otrsub,otrsub]);
          meanInter <- mean(concensusMat[whosub,otrsub]);
          sdInter <- sd(concensusMat[whosub,otrsub]);
          distance <- (meanotr - meanInter)/sqrt((stdotr^2+sdInter^2)/2);
#          cat(c(cid,oid,distance),"\n")
          if (distance < Zthr)
          {
            classID[classID==oid] <- cid;
            sclases <- sclases[!(sclases %in% oid)];
          }
        }
      }
    }
  }
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
  orgnames <-  names(classID);
  classID <- as.numeric(classID);
  names(classID) <- orgnames;
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
  attr(classID,"concensusMat") <- concensusMat
  attr(classID,"pointJaccard") <- pointJaccard
  class(classID) <- "ConsesusLables"
  return (classID);
}


