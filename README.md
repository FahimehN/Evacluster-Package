# Evacluster: Stocastic Algorithms for Cluster Evaluation

[![License](https://img.shields.io/badge/license-GPL--3-blue)](https://www.gnu.org/licenses/gpl-3.0) [![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue)](https://www.r-project.org/) [![CRAN](https://img.shields.io/cran/v/Evacluster)](https://CRAN.R-project.org/package=Evacluster)

Evacluster is an R package that provides stocastic algorithms for evaluating clustering results via consensus clustering and determining optimal cluster numbers. It implements various validity indices and optimization techniques to assess clustering quality.

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
result <- clusterStability(df,clustermethod=MeanShiftCluster)
clusterLabels <- getConsensusCluster(result)
```

## **Available Validity Indices**

The package supports several cluster validity indices including:

-   Silhouette index

-   Dunn index

-   Jacckard.

## **License**

This package is released under the GPL-3 license.
