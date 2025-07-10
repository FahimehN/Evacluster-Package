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
      message(minLables,":",sprintf("%5.3f",lthr),": ",totlabeled,": \n")
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
  attr(classID,"concensusMat") <- concensusMat
  attr(classID,"pointJaccard") <- pointJaccard
  class(classID) <- "ConsesusLables"
  return (classID);
}

#' ConsesusLables summary function
#'
#' This function prints the main statistics of the cluster consensus.
#'
#' @param object A returned object of ConsesusLables function
#' @return the table summarizing the mean and standard deviations of the clustering
#'
#' @export
summary.ConsesusLables<- function(object)
{
  result <- cbind(Mean_intra=attr(object,"Mean_in"),
                 Std_intra=attr(object,"SD_in"),
                 Mean_inter=attr(object,"Mean_out"),
                 Std_inter=attr(object,"SD_out"),
                 Quality=attr(object,"Quality"))
  result <- as.data.frame(result)
  rownames(result) <- names(attr(object,"Mean_in"))
  return (result)
}
  



#' ConsesusLables plot function
#'
#' This function  plots the heat map of the consensus matrix.
#'
#' @param object A returned object of ConsesusLables function
#' @param \dots Additional arguments passed to the heatmap plot function.
#' @return the heatmap and barplot
#'
#' @export
plot.ConsesusLables <- function(LablesResult,...)
{
  op <- par(no.readonly=TRUE)
  
  mycolors <- c("red","green","blue","yellow","orange",
                "red","green","blue","yellow","orange",
                "red","green","blue","yellow","orange",
                "red","green","blue","yellow","orange",
                "red","green","blue","yellow","orange")
  
  
  quality <- attr(LablesResult,"Quality")
  ordermatrix <- attr(LablesResult,"concensusMat")
  theJaccard <- attr(LablesResult,"pointJaccard")
  if (length(LablesResult)>1000)
  {
    LablesResult <- LablesResult[sample(length(LablesResult),1000)]
  }
  ordermatrix <- ordermatrix[names(LablesResult),names(LablesResult)]
  theJaccard <- theJaccard[names(LablesResult)]
  
  orderindex <- 10*LablesResult - theJaccard
  
  orderindex <- order(orderindex)
  ordermatrix <- ordermatrix[orderindex,orderindex]
  rowcolors <- mycolors[clusterLabels]
  rowcolors <- rowcolors[orderindex]
  
  #                             lhei = c(1.5,0.2,4.0),
  
  lmat = rbind(c(6,0,5),c(0,0,2),c(4,1,3))
  lhei = c(2.0,0.2,3.5)
  lwid = c(1.5,0.2,4.0)
  
  hplot <- gplots::heatmap.2(as.matrix(ordermatrix),
                             Rowv=FALSE,Colv=FALSE,
                             RowSideColors = rowcolors,
                             ColSideColors = rowcolors,
                             dendrogram = "none",
                             trace="none",
                             lmat=lmat,
                             lhei=lhei,
                             lwid=lwid,
                             ...)
#  bp <- NULL
  #  print(hplot)
  omd1 <- 0.05+hplot$layout$lwid[1]/sum(hplot$layout$lwid)
  omd1 <- c(omd1,0.90,1.05-hplot$layout$lhei[1]/sum(hplot$layout$lhei))
  omd1 <- c(omd1,0.88)
  par(new=TRUE,cex=0.5,mai=c(0.0,0.40,0.0,0),omd=omd1)
  bp <- barplot(quality,ylab ="Quality",ylim=c(0,1.0))
  par(op)
  
  result <- list(hplot=hplot,bp=bp)
  return(result)
}
  
