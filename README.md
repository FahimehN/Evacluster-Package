# Evacluster: Stochastic Algorithms for Cluster Evaluation

[![License](https://img.shields.io/badge/license-GPL--3-blue)](https://www.gnu.org/licenses/gpl-3.0) [![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue)](https://www.r-project.org/) [![CRAN](https://img.shields.io/cran/v/Evacluster)](https://CRAN.R-project.org/package=Evacluster)

Evacluster is an R package that uses a stochastic algorithms for evaluating clustering results via consensus analysis and determining optimal cluster numbers. It implements various validity indices and to assess clustering quality.

## Installation

You can install the released version from CRAN with:

``` r
install.packages("Evacluster")

# install.packages("devtools")
devtools::install_github("FahimehN/Evacluster-Package")
```

## **Key Features**

-   Stochastic algorithm-based evaluation

-   Allows dimension reduction as a pre-processing step

-   Implements multiple cluster validity indices

-   Supports various distance metrics

-   Visualization tools for co-association evaluation

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
                           trainFraction = 0.95)
clusterLabels <- getConsensusCluster(result)

# Show the clustering stability results
mycolors <- c("red","green","blue","yellow","orange")

table(clusterLabels)

ordermatrix <- result$dataConcensus
 
orderindex <- 10*clusterLabels + result$jaccardpoint
 
orderindex <- order(orderindex)
ordermatrix <- ordermatrix[orderindex,orderindex]
rowcolors <- mycolors[1+clusterLabels]
rowcolors <- rowcolors[orderindex]
 
 
hplot <- gplots::heatmap.2(as.matrix(ordermatrix),
                            Rowv=FALSE,Colv=FALSE,
                            RowSideColors = rowcolors,
                            ColSideColors = rowcolors,
                            dendrogram = "none",
                            trace="none",
                            main="Cluster Co-Association \n (Mean Shift)")
                            
```

## **Available Validity Indices**

The package provides several cluster validity indices including:

-   randIndex

-   Jacckard

-   PAC

## **License**

This package is released under the GPL-3 license.
