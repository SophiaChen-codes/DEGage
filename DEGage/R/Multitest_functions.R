#Alicia Petrany, draft as of 4/25/2023
#This file contains functions for DEGage complete, w
#hich is a full workflow from raw counts to annotated
#groups with differentially expressed genes



#'print.labels
#'
#'prints all labels identified after preprocessing
#'
#'@param labels cell type annotations for every cell
#'@returns nothing
print.labels <- function(labels){
  message("The following cell type annotations were generated:")
  for(i in levels(labels)){
    message(i)
  }
}

#'run.multitest
#'
#'calls DEGage for every comparison necessary in the DEGage_complete workflow
#'
#'@param allcounts the entire matrix of counts for all cells
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param perm.preprocess A logical indicating whether or not to perform the permutation prefiltering step. Setting it to TRUE increases runtimes.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@param perm.pval P value for the permutation test to filter out genes with
#'@param ncores the number of cores to allocated for parallel computing
#'@param writing.dir A directory to write the results of each comparison to. It is HIGHLY reccomended to provide a directory.
#'@return Nothing if a directory is provided. If no directory is provided,
run.multitest <- function(allcounts, group,
                       perm.preprocess,
                       gene.filter.threshold,
                       nperms, nsubsample,
                       perm.pval, ncores, writing.dir){

  if(is.null(writing.dir)){
    outputlist <- list()
  }

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
