#'@title run.multitest
#'@description
#'calls DEGage for every comparison necessary in the DEGage_complete workflow
#'
#'@param allcounts the entire matrix of counts for all cells
#'@param counts A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param perm.preprocess A logical indicating whether or not to perform the permutation prefiltering step. Setting it to TRUE increases runtimes.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@param nsubsample The number of cells to subsample for each condition
#'@param perm.pval P value for the permutation test to filter out genes with
#'@param ncores the number of cores to allocated for parallel computing
#'@param maxiter The maxiumum number of iterations to perform while calculating the cdf
#'@param mean.ratio The minimum ratio between the count means for each group for a gene with an incalculable p-value to be considered differentially expressed
#'@param subsampled.k If true the subsampling procedure is used to estimate k. If false, the random assignment procedure is used.
#'@return Nothing if a directory is provided.
run.multitest <- function(allcounts,
                          group,
                          perm.preprocess,
                          gene.filter.threshold,
                          nperms,
                          nsubsample,
                          perm.pval,
                          ncores,
                          maxiter,
                          mean.ratio,
                          subsampled.k,
                          writing.dir){

  if(is.null(writing.dir)){
    warning("No writing directory was provided. Providing one is strongly recommended")
    outputlist <- list()
  }else{
    if(!dir.exists(writing.dir)){
      stop("Provided writing directory does not exist")
    }
  }

  allcounts <- parse.input(allcounts)

  for(i in 1:(length(levels(group))-1)){
      for(j in (i+1):length(levels(group))){

          message(paste("Comparing", levels(group)[i], "v",
                         levels(group)[j]))

          counts <- cbind(allcounts[,which(group == levels(group)[i])],
                          allcounts[,which(group == levels(group)[j])])
          counts <- as.data.frame(counts)

          tempgroup <- factor(c(rep(0, length(which(group == levels(group)[i]))),
                         rep(1, length(which(group == levels(group)[j])))))

          res <- DEGage(counts, tempgroup, perm.preprocess,
                        gene.filter.threshold, nperms,
                        nsubsample, perm.pval, ncores)

          trialname <- paste(levels(group)[i], "v",
                            levels(group)[j], sep = "")

          if(is.null(writing.dir)){
            outputlist[[(length(outputlist)+1)]] <- res
            names(outputlist)[length(outputlist)] <- trialname
          }
          else{
            write.csv(res, paste(writing.dir, trialname, ".csv", sep = ""))
          }

      }
  }

  if(is.null(writing.dir)){
    return(outputlist)
  }
  else{
    return(NULL)
  }
}
