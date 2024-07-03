#Alicia Petrany, draft as of 5.4.2024
#This file contains main functions
#Install on machine with devtools::install()

#'@title DEGage
#'@description
#'Tests for pairwise differential expression between two count distributions
#'
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
#'@return A dataframe containing regression parameters and p-values for each gene
#'@export
DEGage <- function(counts,
                   group,
                   perm.preprocess = TRUE,
                   gene.filter.threshold = 0.9,
                   nperms = 2000,
                   nsubsample = NA,
                   perm.pval = 0.1,
                   ncores = 4,
                   maxiter = 100,
                   mean.ratio = 4,
                   subsampled.k = T){

  counts <- parse.input(counts)
  #Quality control of user input
  QC(counts, group, gene.filter.threshold ,nperms)

  counts <- counts[rowSums(counts == 0)/ncol(counts) < gene.filter.threshold,] #Filters out genes with all zeros
  outputdf <- data.frame()

  #subsampling for large datasets
  subsampleres <- subsampling(counts, group, nsubsample)
  counts <- subsampleres[[1]]
  group <- subsampleres[[2]]
  rm(subsampleres)

  #configures cores for parallel computing
  cl <- core.config(ncores, nperms, counts, group)

  #calculates permutation test pvalues, then filters counts if perm.preproces == TRUE
  permresultdf <- run.permtest(counts, group, perm.pval, nperms, cl)

  if(perm.preprocess){
   counts <- counts[permresultdf$pval  < perm.pval,]
  }

  #Performs genewise negative binomial regression, generates the df that is to be output with regression parameters
  outputdf <- run.NB.fitting(counts, group, cl)
  outputdf <- run.ZINB.refitting(counts, group, outputdf, cl)

  if(subsampled.k){
    message("calculating K - subsampling")
    outputdf$k <- as.numeric(apply(counts, FUN = calculate_k_subsample, group = group, MARGIN = 1))
  }else{
    message("calculating K - random pairing")
    outputdf$k <- as.numeric(apply(counts, FUN = calculate_k_random, group = group, MARGIN = 1))
  }

  rm(counts)

  #Calculates pvals based off of regression parameters calculated above, adds new column to outputdf containing pvals
  message("Calculating pvals")
  outputdf$pval <- cdf_facilitator(outputdf, maxiter)


  #Adds permutation pvals as output to rows that have not been filtered out
  outputdf$permPvals <- permresultdf$pval[which(permresultdf$gene %in%
                                                rownames(outputdf))]

  stopCluster(cl) #ends parallel computing

  #P value correction for multiple tests
  outputdf$FDR <- p.adjust(outputdf$pval, method = "fdr")

  outputdf <- fix_infs(outputdf, mean.ratio)
  outputdf$pval[is.infinite(outputdf$pval)] <- 1
  outputdf <- organize_output(outputdf)
  return(outputdf)
}

#'@title DEGage_preprocess
#'@description
#'A simple preprocessing workflow to preprocess and generate SingleR cell type
#'annotations for scRNA-seq data.
#'
#'@param counts A dataframe of counts
#'@param min.transcripts.per.gene The minimum number of transcripts a gene can have
#'@param max.transcripts.per.gene The maximum number of transcripts a gene can have
#'@param min.transcripts.per.cell The the minimum number of transcripts required for a cell to be retained
#'@param max.transcripts.per.cell The maximum number of transcripts required for a cell to be retained
#'@param min.cells.per.gene The number of cells a gene must be present in to be retained
#'@param max.dropout.prop Removes genes with zero counts present in a number of cells greater than the specified proportion
#'@param max.mt.percent The maximum proportion of transcripts from mitochondrial genes that can be present in a cell
#'@param mt.prefix The prefix for mitochondrial genes
#'@param max.rRna.precent The maximum proportion of transcripts from ribosomal RNAs that can be present in a cell
#'@param rRna.prefix the prefix for ribosomal genes
#'@return A dataframe of filtered counts
#'@export
DEGage_preprocess <- function(counts,
                              min.transcripts.per.gene = 200,
                              max.transcripts.per.gene = 10000,
                              min.transcripts.per.cell = 500,
                              max.transcripts.per.cell = 20000,
                              min.cells.per.gene = 3,
                              max.dropout.prop = 0.8,
                              max.mt.percent = 0.2,
                              mt.prefix = 'MT-',
                              max.rRna.percent = 0.2,
                              rRna.prefix = 'rr'){

  counts <- counts[rowSums(counts) > min.transcripts.per.gene,]
  counts <- counts[rowSums(counts) < max.transcripts.per.gene,]
  counts <- counts[,colSums(counts) > min.transcripts.per.cell]
  counts <- counts[,colSums(counts) < max.transcripts.per.cell]
  counts <- counts[as.numeric(apply(counts, 1, function(x) sum(x != 0))) > min.cells.per.gene,]
  dropout.props <- as.numeric(apply(counts, 1, function(x) sum(x == 0)))/ncol(counts)
  counts <- counts[dropout.props < max.dropout.prop,]
  mt.percent <- colSums(counts[grepl(mt.prefix, rownames(counts)),])/colSums(counts)
  counts <- counts[,as.numeric(mt.percent) < max.mt.percent]
  rr.percent <- colSums(counts[grepl(rRna.prefix, rownames(counts)),])/colSums(counts)
  counts <- counts[,as.numeric(rr.percent) < max.rRna.percent]
  return(counts)
}
#'@title DEGage_multitest
#'@description
#'Tests for differential expression across more than two conditions
#'
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
#'@param writing.dir A directory to write the results of each comparison to. It is HIGHLY reccomended to provide a directory.
#'@return If a writing directory is not provided, a list of dataframes containing all comparisons is returned
#'@export
DEGage_multitest <- function(counts,
                             group,
                             perm.preprocess = TRUE,
                             gene.filter.threshold = 0.9,
                             nperms = 2000,
                             nsubsample = NA,
                             perm.pval = 0.1,
                             ncores = 4,
                             maxiter = 100,
                             mean.ratio = 4,
                             subsampled.k = T,
                             writing.dir = NULL){

  #add QC's to check for number of group levels
  x <- run.multitest(counts,
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
                     writing.dir = NULL)
  return(x)
}


