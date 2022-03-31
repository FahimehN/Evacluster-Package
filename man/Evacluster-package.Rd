\name{Evacluster-package}
\alias{Evacluster-package}
\alias{Evacluster}
\docType{package}
\title{Evaluation Clustering Methods for Disease Subtypes Diagnosis (Evacluster)}
\description{
Contains a set of clustering methods and evaluation metrics to select the best number of the clusters based on clustering stability.
}

\details{
    \tabular{ll}{
        Package: \tab Evacluster\cr
        Type: \tab Package\cr
        Version: \tab 0.1.0\cr
        Date: \tab 2022-03-25\cr
        License: \tab LGPL (>= 2)\cr
    }
    Purpose: The design of clustering models and evaluation metrics for finding the cluster's number via computing clustering stability.
The best number of clusters is selected via consensus clustering and clustering stability. 
}
\references{Nezhadmoghadam, Fahimeh, and Jose Tamez-Pena. "Risk profiles for negative and positive COVID-19 hospitalized patients.(2021) \emph{Computers in biology and medicine} 136 : 104753.\cr
Fahimeh Nezhadmoghadam, et al., Robust Discovery of Mild Cognitive impairment subtypes and their Risk of Alzheimer's Disease conversion using unsupervised machine learning and Gaussian Mixture Modeling (2021), \emph{Current Alzheimer Research}, 18 (7), 595-606.}
\author{
Fahimeh Nezhadmoghadam, Jose Gerardo Tamez-Pena,
Maintainer: <f.nejad.moghadam@gmail.com>
}
\examples{
    \dontrun{
    ### Evacluster Package Examples ####
    library(datasets)
    data(iris)

   # Split data to training set and testing set
   rndSamples <- sample(nrow(iris),100)
   trainData <- iris[rndSamples,]
   testData <- iris[-rndSamples,]

  
   ## Expectation Maximization Clustering
   # Perform Expectation Maximization Clustering on training set with 3 clusters 
   clsut <- EMCluster(trainData[,1:4],3)
   
   # Predict the labels of the cluster for new data based on cluster labels of the training set
   pre <- predict(clsut,testData[,1:4])
   
   
   ## Fuzzy C-means Clustering
   # Perform Fuzzy C-means Clustering on training set with 3 clusters 
   clsut <- FuzzyCluster(trainData[,1:4],3)
   
   # Predict the labels of the new data 
   pre <- predict(clsut,testData[,1:4])
   
   
   ## hierarchical clustering
   # Perform hierarchical clustering on training set with 3 clusters 
   clsut <- hierarchicalCluster(trainData[,1:4],distmethod="euclidean",clusters=3)
   
   # Predict the labels of the new data 
   pre <- predict(clsut,testData[,1:4])
   
   
   ## K-means Clustering
   # Perform K-means Clustering on training set with 3 clusters 
   clsut <- kmeansCluster(trainData[,1:4],3)
   
   # Predict the labels of the new data 
   pre <- predict(clsut,testData[,1:4])
   
   
   ## Partitioning Around Medoids (PAM) Clustering
   # Perform pam Clustering on training set with 3 clusters 
   clsut <- pamCluster(trainData[,1:4],3)
   
   # Predict the labels of the new data 
   pre <- predict(clsut,testData[,1:4])
   
   
   ## Non-negative matrix factorization (NMF)
   # Perform nmf Clustering on training set with 3 clusters 
   clsut <- nmfCluster(trainData[,1:4],rank=3)
   
   # Predict the labels of the new data 
   pre <- predict(clsut,testData[,1:4])
   
   
   ## t-Distributed Stochastic Neighbor Embedding (t-SNE)
   
   library(mlbench)
   data(Sonar)
 
   rndSamples <- sample(nrow(Sonar),150)
   trainData <- Sonar[rndSamples,]
   testData <- Sonar[-rndSamples,]
 
   # Perform tSNE dimensionality reduction method on training data 
   tsne_trainData <- tsneReductor(trainData[,1:60],dim = 3,perplexity = 10,max_iter = 1000)
   
   # performs an embedding of new data using an existing embedding
   tsne_testData <- predict(tsne_trainData,k=3,testData[,1:60])
   
   
   ## clustering stability function
   # Compute the stability of clustering to select the best number of clusters.
   library(mlbench)
   data(Sonar)
 
   Sonar$Class <- as.numeric(Sonar$Class)
   Sonar$Class[Sonar$Class == 1] <- 0
   Sonar$Class[Sonar$Class == 2] <- 1
   
   # Compute the stability of clustering using kmeans clustering, UMAP as 
   dimensionality reduction method, and feature selection technique
   
  ClustStab <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
                              n_components = 3,featureselection="yes", outcome="Class",
                              fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=3)
   
   
   # Get the labels of the subjects that share the same connectivity
   clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.30,-0.1))


     # Compute the stability of clustering using PAM clustering, tSNE as
     dimensionality reduction method, and feature selection technique
     
   ClustStab <- clusterStability(data=Sonar, clustermethod=pamCluster, dimenreducmethod="tSNE",
                              n_components = 3, perplexity=10,max_iter=100,k_neighbor=2,
                             featureselection="yes", outcome="Class",fs.pvalue = 0.05,
                               randomTests = 100,trainFraction = 0.7,k=3)
          
    # Get the labels of the subjects that share the same connectivity
   clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.30,-0.1))
                     
                     
    # Compute the stability of clustering using hierarchical clustering,
    PCA as dimensionality reduction method, and without applying feature selection
                                 
   ClustStab <- clusterStability(data=Sonar, clustermethod=hierarchicalCluster, 
                               dimenreducmethod="PCA", n_components = 3,featureselection="no",
                               randomTests = 100,trainFraction = 0.7,distmethod="euclidean", 
                               clusters=3)
                               
 # Get the labels of the subjects that share the same connectivity
   clusterLabels <- getConsensusCluster(ClustStab,who="training",thr=seq(0.80,0.30,-0.1))
   
   
   # Show the clustering stability resuldts
   mycolors <- c("red","green","blue","yellow")
 
   ordermatrix <- ClustStab$dataConcensus
 
   heatmapsubsample <- sample(nrow(ordermatrix),70)
 
   orderindex <- 10*clusterLabels + ClustStab$trainJaccardpoint
 
   orderindex <- orderindex[heatmapsubsample]
   orderindex <- order(orderindex)
   ordermatrix <- ordermatrix[heatmapsubsample,heatmapsubsample]
   ordermatrix <- ordermatrix[orderindex,orderindex]
   rowcolors <- mycolors[1+clusterLabels[heatmapsubsample]]
   rowcolors <- rowcolors[orderindex]
 
 
   hplot <- gplots::heatmap.2(as.matrix(ordermatrix),Rowv=FALSE,Colv=FALSE,
                            RowSideColors = rowcolors,ColSideColors = rowcolors,dendrogram = "none",
                            trace="none",main="Cluster Co-Association \n (k=3)")
                            
   
   # Compare the PAC values of clustering stability with different numbers of clusters 
   
   ClustStab2 <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
                              n_components = 3,featureselection="yes", outcome="Class",
                              fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=2)
 
   ClustStab3 <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
                                n_components = 3,featureselection="yes", outcome="Class",
                                fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=3)
 
   ClustStab4 <- clusterStability(data=Sonar, clustermethod=kmeansCluster, dimenreducmethod="UMAP",
                                n_components = 3,featureselection="yes", outcome="Class",
                                fs.pvalue = 0.05,randomTests = 100,trainFraction = 0.7,center=4)
                                
                                
   color_range<- c(black="#FDFC74", orange="#76FF7A", skyblue="#B2EC5D")
 
 
   max.temp <- c(ClustStab2$PAC,ClustStab3$PAC,ClustStab4$PAC) 
 
   barplot(max.temp,xlab = "Number of clusters",ylab = "PAC", names.arg = c( "2","3","4"), 
          ylim=c(0,0.3),col= color_range[1:length(c(1,6,2,6,1))])
                            
   }
}
\keyword{package}
