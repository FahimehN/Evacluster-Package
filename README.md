# Evacluster: Stochastic Algorithms for Cluster Evaluation

[![License](https://img.shields.io/badge/license-GPL--3-blue)](https://www.gnu.org/licenses/gpl-3.0) [![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue)](https://www.r-project.org/) [![CRAN](https://img.shields.io/cran/v/Evacluster)](https://CRAN.R-project.org/package=Evacluster)

Evacluster is an R package that uses a stochastic algorithms for evaluating clustering results via consensus analysis to determining the optimal clustering. It computes various validity indices and to assess clustering quality via the computation of the proportion of ambiguous clustering (PAC) score

## Installation

You can install the released version from CRAN with:

``` r
install.packages("Evacluster")

# install.packages("devtools")
devtools::install_github("FahimehN/Evacluster-Package")
```

## **Key Features**

-   Stochastic-based evaluation

-   Allows dimension reduction as a pre-processing step:

    -   It can be filtered-based: outcome association

    -   It can be based on PCA, tsne, or UMAP

-   Evaluates cluster quality via rand index, Jaccard, and PAC

-   It can be used to get the final Consensus Clusters

## **Basic Usage**

``` r
library(Evacluster)

# Example using iris dataset
data(iris)
df <- iris[, -5]  # Remove species column

# Evaluate clustering with default parameters
result <- clusterStability(df,
                           clustermethod=MeanShiftCluster,
                           featureselection="no",
                           randomTests = 100,
                           trainFraction = 0.5)
#                           nNeighbors=round(0.35*nrow(data)),
#                           bandwidth=rep(0.33,NCOL(df)))

print(result$PAC) # Print the  proportion of ambiguous clustering (PAC)
hist(result$jaccardpoint) # The histogram of the jaccard
hist(result$randIndex)

clusterLabels <- getConsensusCluster(result,thr = seq(0.85, 0.45, -0.1))
barplot(attr(clusterLabels,"Quality"),ylab="Quality",xlab="Cluster",main="Cluster Quality")

table(clusterLabels,iris$Species)


# Show the clustering stability results
mycolors <- c("red","green","blue","yellow","orange")

table(clusterLabels)
plot(iris[,1:2],col=mycolors[clusterLabels],pch=20,cex=2)
text(iris[,1],iris[,2],iris$Species,cex=1)


ordermatrix <- result$dataConcensus
 
orderindex <- 10*clusterLabels + result$jaccardpoint
 
orderindex <- order(orderindex)
ordermatrix <- ordermatrix[orderindex,orderindex]
rowcolors <- mycolors[clusterLabels]
rowcolors <- rowcolors[orderindex]
 
 
hplot <- gplots::heatmap.2(as.matrix(ordermatrix),
                            Rowv=FALSE,Colv=FALSE,
                            RowSideColors = rowcolors,
                            ColSideColors = rowcolors,
                            dendrogram = "none",
                            trace="none",
                            main="Cluster Co-Association \n (Mean Shift)")


## Random numbers should not create a clear co-association matrix

df <- as.data.frame(matrix(rnorm(200 * 3), nrow = 200, ncol = 3))
result <- clusterStability(df,
                           clustermethod=kmeansCluster,
                           featureselection="no",
                           randomTests = 100,
                           trainFraction = 0.5,
                           center=2)
print(result$PAC) # Print the  proportion of ambiguous clustering (PAC)
hist(result$jaccardpoint) # The histogram of the jaccard
hist(result$randIndex)

clusterLabels <- getConsensusCluster(result,who = "testing",thr = seq(0.9, 0.5, -0.1))
barplot(attr(clusterLabels,"Quality"),ylab="Quality",xlab="Cluster",main="Cluster Quality")
clusterLabels <- getConsensusCluster(result,who = "training",thr = seq(0.9, 0.5, -0.1))
barplot(attr(clusterLabels,"Quality"),ylab="Quality",xlab="Cluster",main="Cluster Quality")


table(clusterLabels)
orderindex <- 10*clusterLabels + result$jaccardpoint
ordermatrix <- result$testConsesus
 
orderindex <- order(orderindex)
ordermatrix <- ordermatrix[orderindex,orderindex]
rowcolors <- mycolors[clusterLabels]
rowcolors <- rowcolors[orderindex]
 
 
hplot <- gplots::heatmap.2(as.matrix(ordermatrix),
                            Rowv=FALSE,Colv=FALSE,
                            RowSideColors = rowcolors,
                            ColSideColors = rowcolors,
                            dendrogram = "none",
                            trace="none",
                            main="Cluster Co-Association \n (kmeans)")
                            
```

## **Available Validity Indices**

The package provides several cluster validity indices including:

-   Rand Index

-   Jaccard similarity

-   PAC

## **License**

This package is released under the GPL-3 license.
