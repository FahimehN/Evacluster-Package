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
        Date: \tab 2022-03-24\cr
        License: \tab LGPL (>= 2)\cr
    }
    Purpose: The design of clustering models and evaluation metrics for finding the cluster's number via computing clustering stability.
The best number of clusters is selected via consensus clustering and clustering stability. 
}
\references{Nezhadmoghadam, Fahimeh, and Jose Tamez-Peña. "Risk profiles for negative and positive COVID-19 hospitalized patients.(2021) \emph{Computers in biology and medicine} 136 : 104753.
Fahimeh Nezhadmoghadam, et al., Robust Discovery of Mild Cognitive impairment subtypes and their Risk of Alzheimer's Disease conversion using unsupervised machine learning and Gaussian Mixture Modeling (2021), \emph {Current Alzheimer Research}, 18 (7), 595-606. 
}
\author{
Fahimeh Nezhadmoghadam, Jose Gerardo Tamez-Pena 
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

  
   # Perform Expectation Maximization Clustering of training set with 3 clusters 
   clsut <- EMCluster(trainData[,1:4],3)
   
   # Predict the labels of the cluster for new data based on cluster labels of the training set
   pre <- predict(clsut,testData[,1:4])
   
   
   # Perform Fuzzy C-means Clustering of training set with 3 clusters 
   clsut <- FuzzyCluster(trainData[,1:4],3)
   
   # Predict the labels of the cluster for new data based on cluster labels of the training set
   pre <- predict(clsut,testData[,1:4])
 

   


   

   }
}
\keyword{package}
