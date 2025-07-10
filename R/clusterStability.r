#' clustering stability function
#'
#' This function computes the stability of clustering that helps to select the best number of clusters.
#' Feature selection and dimensionality reduction methods can be used before clustering 
#' the data.
#'
#' @param data A Data set
#' @param clustermethod The clustering method. This can be one of "Mclust","pamCluster","kmeansCluster", "hierarchicalCluster",and "FuzzyCluster".
#' @param dimenreducmethod The dimensionality reduction method. This must be one of "Auto", "none","UMAP","tSNE", and "PCA".
#' @param n_components The dimension of the space that data embed into. It can be set to any integer value in the range of 2 to 100.
#' @param perplexity The Perplexity parameter that determines the optimal number of neighbors in tSNE method.(it is only used in the tSNE reduction method)
#' @param max_iter The maximum number of iterations for performing tSNE reduction method.
#' @param k_neighbor The k_neighbor is used for computing the means of #neighbors with min distance (#Neighbor=sqrt(#Samples/k) for performing an embedding of new data using an existing embedding in the tSNE method.
#' @param featureselection This parameter determines whether feature selection is applied before clustering data or not. if used, it should be "yes", otherwisw "no".
#' @param outcome The outcome feature is used for feature selection.
#' @param fs.pvalue The threshold pvalue used for feature selection process. The default value is 0.05.
#' @param randomTests The number of iterations of the clustering process for computing the cluster stability.
#' @param trainFraction This parameter determines the ratio of training data. The default value is 0.5.
#' @param pac.thr  The pac.thr is the threshold to use for computing the proportion of ambiguous clustering (PAC) score. It is as the fraction of sample pairs with consensus indices falling in the interval.The default value is 0.1.
#' @param plotClustering  if TRUE the class-labeled scatter plot of the first two dimensions will be shown.
#' @param \dots Additional arguments passed to the clustering algorithm.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item randIndex - A vector of the Rand Index that computes a similarity measure between two clusterings. 
#'   \item jaccIndex - A vector of jaccard Index that measures how frequently pairs of items are joined together in two clustering data sets.
#'   \item randomSamples - A vector with indexes of selected samples for training in each iteration.
#'   \item clusterLabels - A vector with clusters' labels in all iterations. jaccardpoint
#'   \item jaccardpoint - The corresponding Jaccard index for each data point of testing set
#'   \item averageNumberofClusters - The mean Number of Clusters.
#'   \item testConsesus - A vector of consensus clustering results of testing set.
#'   \item trainRandIndex - A vector of the Rand Index for training set.
#'   \item trainJaccIndex - A vector of the jaccard Index for training set.
#'   \item trainJaccardpoint - The corresponding Jaccard index for each data point of training set.
#'   \item PAC - The proportion of ambiguous clustering (PAC) score.
#'   \item dataConcensus - A vector of consensus clustering results of training set.
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
#'
#' ClustStab <- clusterStability(data=Sonar, clustermethod=pamCluster, dimenreducmethod="tSNE",
#'                               n_components = 3, perplexity=10,max_iter=100,k_neighbor=2,
#'                               featureselection="yes", outcome="Class",fs.pvalue = 0.05,
#'                               randomTests = 100,trainFraction = 0.7,k=3)
#'
#'
#' ClustStab <- clusterStability(data=Sonar, clustermethod=hierarchicalCluster, 
#'                               dimenreducmethod="PCA", n_components = 3,featureselection="no",
#'                               randomTests = 100,trainFraction = 0.7,distmethod="euclidean",
#'                               clusters=3)
#'
#'}
#' @export
clusterStability <- function(data=NULL, 
                             clustermethod=MeanShiftCluster, 
                             dimenreducmethod=c("Auto","none","PCA","tsne","UMAP"),
                             n_components = 3,
                             perplexity = 25,
                             max_iter = 1000,
                             k_neighbor=3,
                             featureselection=NULL,
                             outcome=NULL,
                             fs.pvalue = 0.05,
                             randomTests = 20,
                             trainFraction = 0.5,
                             pac.thr=0.1,
                             plotClustering=FALSE,
                             ...)
{
  clusterLabels <- list();
  randomSamples <- list();
  numberofClusters <- 0;
  testCounts <- numeric(nrow(data))
  randomSeeds <- sample(randomTests);

  dimenreducmethod <- match.arg(dimenreducmethod);
  if (dimenreducmethod=="Auto")
  {
    if (ncol(data)>5)
    {
      dimenreducmethod <- "UMAP"
    }
    else
    {
      message("No dimension reduction: Using all features\n")
      dimenreducmethod <- NULL
    }
  }
  
  
  for (i in 1:randomTests)
  {
    randomSamples[[i]] <- sample(nrow(data),trainFraction*nrow(data));
    message(paste('iteration= ',i,':',randomTests))
    
    data <- data[, ! names(data) %in% outcome, drop = F]
    tempdata <- data
    
    ### Feature Selection ###
    
    if (!is.null(featureselection))
    {
      if (featureselection == "yes")
      {
        cat("With Feature Selection!!!...")
        message(paste('data Before FS=',nrow(data)))
        
        FS <- names(FRESA.CAD::univariate_Wilcoxon(data = data[randomSamples[[i]],],
                                                   Outcome = outcome,
                                                   pvalue = fs.pvalue))
        tempdata <- data.frame(data[,FS])
        
        message(paste('Number of features= ',length(FS)))
        message("Feature selection was Done!\n")
      }
    }

    ### Dimension Reduction ###
    
    if (!is.null(dimenreducmethod))
    {
      if(dimenreducmethod == "UMAP")
      {
        umapData <- uwot::umap(tempdata[randomSamples[[i]],], ret_model = TRUE,
                               n_components = n_components)
        
        umaptestData <- uwot::umap_transform(data[-randomSamples[[i]],], umapData)
        tempdata <- as.data.frame(matrix(0,nrow(data),ncol(umapData$embedding)));
        rownames(tempdata) <- rownames(data);
        colnames(tempdata) <- colnames(umapData$embedding);
        #        print(c(nrow(tempdata),ncol(tempdata)))
        tempdata[randomSamples[[i]],] <- as.data.frame(umapData$embedding)
        tempdata[-randomSamples[[i]],] <- as.data.frame(umaptestData)
        #        print(c(nrow(tempdata),ncol(tempdata)))
        message("UMAP was Done!")
      }
      else if (dimenreducmethod == "tSNE")
      {
        tsneData <- tsneReductor(tempdata[randomSamples[[i]],],
                                 dim=n_components,perplexity=perplexity,
                                 max_iter=max_iter)
        
        dupIndex <- which(duplicated(tempdata[randomSamples[[i]],], fromLast = TRUE) %in% TRUE)
        
        if(length(dupIndex) != 0)
        {
          randomSamples[[i]] <- randomSamples[[i]][-dupIndex]
        }
        
        tsnetestData <- predict(tsneData,k=k_neighbor,tempdata[-randomSamples[[i]],])
        
        tempdata[randomSamples[[i]],] <- as.data.frame(tsneData$tsneY) 
        tempdata[-randomSamples[[i]],] <- as.data.frame(tsnetestData$tsneY)
        message("t-SNE was Done!")
      }
      else if (dimenreducmethod == "PCA")
      {
        pcaData <- stats::prcomp(tempdata[randomSamples[[i]],],rank. = n_components)
        pcatestData <- stats::predict(pcaData,tempdata[-randomSamples[[i]],])
        tempdata[randomSamples[[i]],] <- as.data.frame(pcaData$x)
        tempdata[-randomSamples[[i]],] <- as.data.frame(pcatestData)
        message("PCA was Done!")
      }
    }
    #    print(c(nrow(tempdata),ncol(tempdata)))
    
    mod1 <- clustermethod(tempdata[randomSamples[[i]],],...);
    clusterLabels[[i]] <- predict(mod1,tempdata); 
    names(clusterLabels[[i]]$classification) <- rownames(tempdata) #data
    collab <- clusterLabels[[i]]$classification[randomSamples[[i]]];
    
    if (plotClustering)
    {
      graphics::plot(tempdata[randomSamples[[i]],1:2],col = collab,main=sprintf("%d",i));
    }
    
    numberofClusters <- numberofClusters + length(table(clusterLabels[[i]]$classification))
    testCounts[-randomSamples[[i]]] <- testCounts[-randomSamples[[i]]] + 1;
    set.seed(randomSeeds[i]);
  }
  
  numberofClusters <- numberofClusters/randomTests;
  message("Done Testing:")
  randIndex <- numeric();
  jaccIndex <- numeric();
  meanJaccard <- numeric();
  jaccardpoint <- numeric(nrow(data));
  jaccardpointcount <- numeric(nrow(data));
  trainrandIndex <- numeric();
  trainjaccIndex <- numeric();
  trainmeanJaccard <- numeric();
  trainjaccardpoint <- numeric(nrow(data));
  trainjaccardpointcount <- numeric(nrow(data));
  for (i in 1:(randomTests - 1))
  {
    for (j in (i + 1):randomTests)
    {
      outsamples <- unique(c(randomSamples[[i]],randomSamples[[j]]))
      if ((nrow(data) - length(outsamples)) > 10)
      {
        randIndex <- c(randIndex,mclust::adjustedRandIndex(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]));
        jaccard <- FRESA.CAD::jaccardMatrix(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]);
        jaccIndex <- c(jaccIndex,jaccard$balancedMeanJaccard);
        meanJaccard <- c(meanJaccard,mean(jaccard$elementJaccard));
        jaccardpoint[-outsamples] <- jaccardpoint[-outsamples] + jaccard$elementJaccard;
        jaccardpointcount[-outsamples] <- jaccardpointcount[-outsamples] + 1;
      }
      insamples <- randomSamples[[i]][randomSamples[[i]] %in% randomSamples[[j]]]
      if ((nrow(data) - length(insamples)) > 10)
      {
        trainrandIndex <- c(trainrandIndex,mclust::adjustedRandIndex(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]));
        trainjaccard <- FRESA.CAD::jaccardMatrix(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]);
        trainjaccIndex <- c(trainjaccIndex,trainjaccard$balancedMeanJaccard);
        trainmeanJaccard <- c(trainmeanJaccard,mean(trainjaccard$elementJaccard));
        trainjaccardpoint[insamples] <- trainjaccardpoint[insamples] + trainjaccard$elementJaccard;
        trainjaccardpointcount[insamples] <- trainjaccardpointcount[insamples] + 1;
      }
    }
  }
  message("After Jacckard:")
  jaccardpoint[jaccardpointcount > 0] <- jaccardpoint[jaccardpointcount > 0]/jaccardpointcount[jaccardpointcount > 0];
  names(jaccardpoint) <- rownames(data);
  trainjaccardpoint[trainjaccardpointcount > 0] <- trainjaccardpoint[trainjaccardpointcount > 0]/trainjaccardpointcount[trainjaccardpointcount > 0];
  names(trainjaccardpoint) <- rownames(data);
  
  testConsesus <- matrix(0,nrow = nrow(data), ncol = nrow(data))
  colnames(testConsesus) <- rownames(data)
  rownames(testConsesus) <- rownames(data)
  countMat <- testConsesus;
  dataConcensus <- testConsesus
  totwts <- 0;
  for (i in 1:randomTests)
  {
    testset <- rownames(data[-randomSamples[[i]],])
    aclassLabels <- clusterLabels[[i]]$classification;
    nclus <- length(table(aclassLabels))
    wts <- (1.0-0.999*(nclus < 2))/(1.0+abs(nclus-numberofClusters));
    classLabels <- aclassLabels[testset];
    btestset <- rownames(data) %in% testset;
    for (id in testset)
    {
      testConsesus[id,btestset] <- testConsesus[id,btestset] + wts*(classLabels == aclassLabels[id]);
      countMat[id,btestset] <- countMat[id,btestset] + wts;
    }
    classLabels <- clusterLabels[[i]]$classification;
    for (id in 1:nrow(data))
    {
      dataConcensus[id,] <- dataConcensus[id,] + wts*(classLabels == classLabels[id]);
    }
    totwts <- totwts + wts;
  }
  message("After Counting.")
  testConsesus[countMat > 0] <- testConsesus[countMat > 0]/countMat[countMat > 0];
  dataConcensus <- dataConcensus/totwts;
  totmatpts <- nrow(data)^2;
  pac <- sum(testConsesus[(testConsesus > pac.thr) & (testConsesus < (1.0 - pac.thr))])/totmatpts;
  pac <- c(pac,sum(testConsesus[(testConsesus > 0.9*pac.thr) & (testConsesus < (1.0 - 0.9*pac.thr))])/totmatpts);
  pac <- c(pac,sum(testConsesus[(testConsesus > 1.1*pac.thr) & (testConsesus < (1.0 - 1.1*pac.thr))])/totmatpts);
  pac <- c(pac,sum(testConsesus[(testConsesus > 0.8*pac.thr) & (testConsesus < (1.0 - 0.8*pac.thr))])/totmatpts);
  pac <- c(pac,sum(testConsesus[(testConsesus > 1.2*pac.thr) & (testConsesus < (1.0 - 1.2*pac.thr))])/totmatpts);
  
  
  result <- list(randIndex = randIndex,jaccIndex = jaccIndex,randomSamples = randomSamples,
                 clusterLabels=clusterLabels,jaccardpoint=jaccardpoint, averageNumberofClusters=numberofClusters,
                 testConsesus=testConsesus,trainRandIndex = trainrandIndex,trainJaccIndex = trainjaccIndex,
                 trainJaccardpoint=trainjaccardpoint,PAC=pac,dataConcensus=dataConcensus);
  class(result) <- "ClusterStability"
  return(result);
}

#' ClusterStability summary function
#'
#' This function prints the main features of the cluster evaluation.
#'
#' @param object A returned object of ClusterStability function
#' @return the table summarizing the mean and standard deviations of the clustering
#'
#' @export
summary.ClusterStability <- function(clusResult)
{
  dc <- clusResult$clusterLabels
  attributes(dc) <- NULL
  numclusters <- unlist(lapply(lapply(dc,table),length))
  AverageClusters <- mean(numclusters)
  stdClusters <- sd(numclusters)
  meanRandTrain <- mean(clusResult$trainRandIndex)
  sdRandTrain <- sd(clusResult$trainRandIndex)
  meanRandTest <- mean(clusResult$randIndex)
  sdRandTest <- sd(clusResult$randIndex)
  meanJaccard <- mean(clusResult$jaccIndex)
  sdJaccard <- sd(clusResult$jaccIndex)
  meanJaccardTrain <- mean(clusResult$trainJaccIndex)
  sdJaccardTrain <- sd(clusResult$trainJaccIndex)
  meanAtPointJaccard <- mean(clusResult$jaccardpoint)
  sdAtPointJaccard <- sd(clusResult$jaccardpoint)
  meanAtPointJaccardTrain <- mean(clusResult$trainJaccardpoint)
  sdAtPointJaccardTrain <- sd(clusResult$trainJaccardpoint)
  meanPAC <- mean(clusResult$PAC)
  sdPAC <- sd(clusResult$PAC)
  
  result <- c(AverageClusters,stdClusters)
  result <- rbind(result,c(meanRandTrain,sdRandTrain))
  result <- rbind(result,c(meanRandTest,sdRandTest))
  result <- rbind(result,c(meanJaccardTrain,sdJaccardTrain))
  result <- rbind(result,c(meanJaccard,sdJaccard))
  result <- rbind(result,c(meanAtPointJaccardTrain,sdAtPointJaccardTrain))
  result <- rbind(result,c(meanAtPointJaccard,sdAtPointJaccard))
  result <- rbind(result,c(meanPAC,sdPAC))
  result <- as.data.frame(result)
#  print(result)
  colnames(result) <- c("Mean","SD")
  rownames(result) <- c("No. Clusters",
                        "Train Rand",
                        "Test Rand",
                        "Train Jaccard",
                        "Test Jaccard",
                        "Train At Point",
                        "Test At Point",
                        "PAC")
  
  return(result)
}


#' ClusterStability plot function
#'
#' This function box plots the indexes of the cluster evaluation.
#'
#' @param object A returned object of ClusterStability function
#' @param \dots Additional arguments passed to the boxplot function.
#' @return the box plots data
#'
#' @export
plot.ClusterStability <- function(clusResult,...)
{
  op <- par(no.readonly=TRUE)
  
  databoxes <- list(Train_Rand=clusResult$trainRandIndex,
                    Rand=clusResult$randIndex,
                    Train_Jaccard=clusResult$trainJaccIndex,
                    Jaccard=clusResult$jaccIndex,
                    Train_At_Point=clusResult$trainJaccardpoint,
                    At_Point=clusResult$jaccardpoint
                    ) 
  par(cex=0.95,mar=c(8,4,4,2)+0.1)
  bp <- boxplot(databoxes,notch=TRUE,las=2,ylab="Score",...)
  par(op)
  return(bp)
}


