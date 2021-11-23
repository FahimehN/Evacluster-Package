clusterOutcomeStability <- function(clusterLabels=NULL,randomSamples=NULL,outcomeLabels=NULL)
{
  outcomRandIndex <- numeric();
  outcomejaccIndex <- numeric();
  outcomemeanJaccard <- numeric();
  for (i in 1:length(clusterLabels))
  {
    if (!is.null(outcomeLabels))
    {
      outcomRandIndex <- c(outcomRandIndex,adjustedRandIndex(clusterLabels[[i]]$classification[-randomSamples[[i]]],outcomeLabels[-randomSamples[[i]]]));
      outcomeJaccard <- jaccardMatrix(clusterLabels[[i]]$classification[-randomSamples[[i]]],outcomeLabels[-randomSamples[[i]]]);
      outcomejaccIndex <- c(outcomejaccIndex,outcomeJaccard$balancedMeanJaccard);
      outcomemeanJaccard <- c(outcomemeanJaccard,mean(outcomeJaccard$elementJaccard));
    }
  }
  
  result <- list(outcomRandIndex=outcomRandIndex,outcomejaccIndex=outcomejaccIndex,outcomemeanJaccard=outcomemeanJaccard);
  class(result) <- "ClusterOutcomeStability"
  return(result);
}
