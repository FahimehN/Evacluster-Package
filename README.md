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
result <- clusterStability(df)

print(summary(result))
pt <- plot(result,main="Evaluation Indexes")

clusterLabels <- getConsensusCluster(result)
print(table(clusterLabels,iris$Species))
print(summary(clusterLabels))
pt <- plot(clusterLabels,main="Consensus Matrix",ylab="Sample",xlab="Sample")

                            
```

## **Available Validity Indices**

The package provides several cluster validity indices including:

-   Rand Index

-   Jaccard similarity

-   PAC

## **License**

This package is released under the GPL-3 license.
