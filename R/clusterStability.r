#' clustering stability function
#'
#' This function computes the stability of clustering that helps to select the best number of clusters.
#' Feature selection and dimensionality reduction methods can be used before clustering 
#' the data.
#'
#' @param data A Data set
#' @param clustermethod The clustering method. This can be one of "Mclust","pamCluster","kmeansCluster", "hierarchicalCluster",and "FuzzyCluster".
#' @param dimenreducmethod The dimensionality reduction method. This must be one of "UMAP","tSNE", and "PCA".
#' @param n_components The dimension of the space that data embed into. It can be set to any integer value in the range of 2 to 100.
#' @param perplexity The Perplexity parameter that determines the optimal number of neighbors in tSNE method.(it is only used in the tSNE reduction method)
#' @param max_iter The maximum number of iterations for performing tSNE reduction method.
#' @param k_neighbor The k_neighbor is used for computing the means of #neighbors with min distance (#Neighbor=sqrt(#Samples/k) for performing an embedding of new data using an existing embedding in the tSNE method.
#' @param featureselection This parameter determines whether feature selection is applied before clustering data or not. if used, it should be "yes", otherwisw "no".
#' @param outcome The outcome feature is used for feature selection.
#' @param fs.pvalue The threshold pvalue used for feature selection process. The default value is 0.05.
#' @param randomTests The number of iterations of the clustering process for computing the cluster stability.
#' @param trainFraction This parameter determines the ratio of training data. The default value is 0.5.
#' @param pac.thr  The pac.thr is the thresold to use for computing the proportion of ambiguous clustering (PAC) score. It is as the fraction of sample pairs with consensus indices falling in the interval.The default value is 0.1.
#' @param \dots Additional arguments passed to clusterStability().
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
clusterStability <- function(data=NULL, clustermethod=NULL, dimenreducmethod=NULL,
                             n_components = 3,perplexity = 25,max_iter = 1000,k_neighbor=3,
                             featureselection=NULL ,outcome=NULL,fs.pvalue = 0.05,
                             randomTests = 20, trainFraction = 0.5,pac.thr=0.1,...)
{
  clusterLabels <- list();
  randomSamples <- list();
  numberofClusters <- 0;
  testCounts <- numeric(nrow(data))
  randomSeeds <- sample(randomTests);
  
  for (i in 1:randomTests)
  {
    randomSamples[[i]] <- sample(nrow(data),trainFraction*nrow(data));
    message(paste('iteration= ',i,':',randomTests))
    
    ### Feature Selection ###
    if (featureselection == "yes")
    {
      message(paste('data Before FS=',nrow(data)))
      
      FS <- names(FRESA.CAD::univariate_Wilcoxon(data = data[randomSamples[[i]],],
                                      Outcome = outcome,
                                      pvalue = fs.pvalue))
      tempdata <- data.frame(data[,FS])
      
      message(paste('Number of features= ',length(FS)))
      warning("Feature selection was Done!")
    }
    else {
      data <- data[, ! names(data) %in% outcome, drop = F]
      tempdata <- data
      warning("Without Feature Selection!!!")
      }
    
    ### Dimension Reduction ###
    
    if (!is.null(dimenreducmethod))
    {
      if(dimenreducmethod == "UMAP")
      {
        umapData <- uwot::umap(tempdata[randomSamples[[i]],], ret_model = TRUE,
                         n_components = n_components)

        umaptestData <- uwot::umap_transform(tempdata[-randomSamples[[i]],], umapData)
        
        tempdata[randomSamples[[i]],] <- as.data.frame(umapData$embedding)
        tempdata[-randomSamples[[i]],] <- as.data.frame(umaptestData)
        warning("UMAP was Done!")
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
        warning("t-SNE was Done!")
      }
      else if (dimenreducmethod == "PCA")
      {
        pcaData <- stats::prcomp(tempdata[randomSamples[[i]],],rank. = n_components)
        pcatestData <- stats::predict(pcaData,tempdata[-randomSamples[[i]],])
        tempdata[randomSamples[[i]],] <- as.data.frame(pcaData$x)
        tempdata[-randomSamples[[i]],] <- as.data.frame(pcatestData)
        warning("PCA was Done!")
      }
      else {warning("Package does not support the selected reduction method !!")}
    }
    
    mod1 <- clustermethod(tempdata[randomSamples[[i]],],...);
    clusterLabels[[i]] <- predict(mod1,tempdata); 
    names(clusterLabels[[i]]$classification) <- rownames(tempdata) #data
    graphics::plot(data[,1:2],col = clusterLabels[[i]]$classification,main=sprintf("%d",i));
    numberofClusters <- numberofClusters + length(table(clusterLabels[[i]]$classification))
    testCounts[-randomSamples[[i]]] <- testCounts[-randomSamples[[i]]] + 1;
    set.seed(randomSeeds[i]);
  }

  numberofClusters <- numberofClusters/randomTests;
  warning("Done Testing:")
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
  warning("After Jacckard:")
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
    wts <- (1.0-0.99*(nclus < 2))/(1.0+abs(nclus-numberofClusters));
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
  warning("After Counting.")
  testConsesus[countMat > 0] <- testConsesus[countMat > 0]/countMat[countMat > 0];
  dataConcensus <- dataConcensus/totwts;
  pac <- sum(testConsesus[(testConsesus > pac.thr) & (testConsesus < (1.0 - pac.thr))])/nrow(data)/nrow(data);


  result <- list(randIndex = randIndex,jaccIndex = jaccIndex,randomSamples = randomSamples,
                 clusterLabels=clusterLabels,jaccardpoint=jaccardpoint, averageNumberofClusters=numberofClusters,
                 testConsesus=testConsesus,trainRandIndex = trainrandIndex,trainJaccIndex = trainjaccIndex,
                 trainJaccardpoint=trainjaccardpoint,PAC=pac,dataConcensus=dataConcensus);
  class(result) <- "ClusterStability"
  return(result);
}


