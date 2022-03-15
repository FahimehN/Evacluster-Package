#' clustering stability function
#'
#' This function computes the clustering stability to select the best number of clusters.
#'  Feature selection and dimensionality reduction methods can be used before clustering 
#'  the data.
#'
#' @param data A Data set
#' @param k The number of Clusters
#' @return A list of cluster labels and a R object of class "fcm {ppclust}"
#' @examples
#' library(datasets)
#' data(iris)
#'
#' rndSamples <- sample(nrow(iris),100)
#' trainData <- iris[rndSamples,]
#' testData <- iris[-rndSamples,]
#'
#' cls <- FuzzyCluster(trainData[,1:4],3)
#' @export
clusterStability <- function(data=NULL, clustermethod=NULL, dimenreducmethod=NULL,
                             n_components = 3,perplexity = 25,max_iter = 1000,k_neighbor=3,
                             featureselection=NULL ,outcome=NULL,fs.pvalue = 0.05, 
                             randomTests = 20, trainFraction = 0.5,pac.thr=0.1, ...)
{
  clusterLabels <- list();
  randomSamples <- list();
  numberofClusters <- 0;
  testCounts <- numeric(nrow(data))
  randomSeeds <- sample(randomTests);
  
  for (i in 1:randomTests)
  {
    randomSamples[[i]] <- sample(nrow(data),trainFraction*nrow(data));
    message(paste('Done= ',i))
    
    ### Feature Selection ###
    if (featureselection == "YES")
    {
      message(paste('data Before FS=',nrow(data)))
      
      FS <- names(univariate_Wilcoxon(data = data[randomSamples[[i]],],
                                      Outcome = outcome,
                                      pvalue = fs.pvalue))
      tempdata <- data.frame(data[,FS])
      
      message(paste('Number of features= ',length(FS)))
      cat("Feature selection was Done!")
    }
    else {
      data <- data[, ! names(data) %in% outcome, drop = F]
      tempdata <- data
      cat("Without Feature Selection!!!")
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
        cat("UMAP was Done!")
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
        cat("t-SNE was Done!")
      }
      else if (dimenreducmethod == "PCA")
      {
        pcaData <- stats::prcomp(tempdata[randomSamples[[i]],],rank. = n_components)
        pcatestData <- predict(pcaData,tempdata[-randomSamples[[i]],])
        tempdata[randomSamples[[i]],] <- as.data.frame(pcaData$x)
        tempdata[-randomSamples[[i]],] <- as.data.frame(pcatestData)
        cat("PCA was Done!")
      }
      else {cat("Package does not support the selected reduction method !!")}
    }
    
    mod1 <- clustermethod(tempdata[randomSamples[[i]],],...);
    clusterLabels[[i]] <- predict(mod1,tempdata); #data[-randomSamples[[i]] #tempdata[-randomSamples[[i]],]
    names(clusterLabels[[i]]$classification) <- rownames(tempdata) #data
    plot(data[,1:2],col = clusterLabels[[i]]$classification,main=sprintf("%d",i));
    numberofClusters <- numberofClusters + length(table(clusterLabels[[i]]$classification))
    testCounts[-randomSamples[[i]]] <- testCounts[-randomSamples[[i]]] + 1;
    set.seed(randomSeeds[i]);
  }

  numberofClusters <- numberofClusters/randomTests;
  cat("Done Testing:")
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
        randIndex <- c(randIndex,adjustedRandIndex(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]));
        jaccard <- jaccardMatrix(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]);
        jaccIndex <- c(jaccIndex,jaccard$balancedMeanJaccard);
        meanJaccard <- c(meanJaccard,mean(jaccard$elementJaccard));
        jaccardpoint[-outsamples] <- jaccardpoint[-outsamples] + jaccard$elementJaccard;
        jaccardpointcount[-outsamples] <- jaccardpointcount[-outsamples] + 1;
      }
      insamples <- randomSamples[[i]][randomSamples[[i]] %in% randomSamples[[j]]]
      if ((nrow(data) - length(insamples)) > 10)
      {
        trainrandIndex <- c(trainrandIndex,adjustedRandIndex(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]));
        trainjaccard <- jaccardMatrix(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]);
        trainjaccIndex <- c(trainjaccIndex,trainjaccard$balancedMeanJaccard);
        trainmeanJaccard <- c(trainmeanJaccard,mean(trainjaccard$elementJaccard));
        trainjaccardpoint[insamples] <- trainjaccardpoint[insamples] + trainjaccard$elementJaccard;
        trainjaccardpointcount[insamples] <- trainjaccardpointcount[insamples] + 1;
      }
    }
  }
  cat("After Jacckard:")
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
  cat("After Counting.")
  testConsesus[countMat > 0] <- testConsesus[countMat > 0]/countMat[countMat > 0];
  dataConcensus <- dataConcensus/totwts;
  pac <- sum(testConsesus[(testConsesus > pac.thr) & (testConsesus < (1.0 - pac.thr))])/nrow(data)/nrow(data);


  result <- list(randIndex = randIndex,jaccIndex = jaccIndex,meanJaccard = meanJaccard,randomSamples = randomSamples,
                 clusterLabels=clusterLabels, jaccardpoint=jaccardpoint,averageNumberofClusters=numberofClusters,
                 testConsesus=testConsesus,trainRandIndex = trainrandIndex,trainJaccIndex = trainjaccIndex,trainMeanJaccard = trainmeanJaccard,
                 trainJaccardpoint=trainjaccardpoint,PAC=pac,dataConcensus=dataConcensus);
  class(result) <- "ClusterStability"
  return(result);
}


plot.ClusterStability <- function(object,...)
{
  plot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
  boxplot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}

summary.ClusterStability <- function(object,...)
{
  summary(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}
