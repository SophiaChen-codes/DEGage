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
DEGage <- function(counts, group, perm.preprocess = TRUE,
                   gene.filter.threshold = 1, nperms = 2000,
                   nsubsample = NA, perm.pval = 0.1, ncores = 4,
                   maxiter = 100, mean.ratio = 1.4, subsampled.k = T){

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
#'@param input either a dataframe or matrix containing counts, or a path to an mtx directory
#'@param dir.type the type of directory that a path, if input, is directing to. Can be 'mtx'.
#'@param min.nFeatureRNA The minimum number of genes that must be present in a cell
#'@param max.nFeatureRNA The maxiumum number of genes that can be present in a cell
#'@param mt.percent A value between 0 and 1 indicating the maximum proporition of mitochondrial genes allowed in a cell
#'@param cell.annotations A boolean indicating whether or not to perform celltype annotations.
#'@param celltype.min the minimum number of cells of a given cell type to count as valid.
#'@param output Format of the output. Can be "Seurat" or "df"
#'@return A data structure of designated user input that contains count information
#'@export
DEGage_preprocess <- function(input, dir.type = 'mtx', min.nFeatureRNA = 200,
                              max.nFeatureRNA = 8000, mt.percent = .2,
                              cell.annotations = TRUE,
                              celltype.min = 20, output = "Seurat"){

  preprocess.QC(input, dir.type, output ='df')

  sobj <- parse.input.preprocess(input, dir.type)
  sobj <- filter.and.normalize(sobj, min.nFeatureRNA, max.nFeatureRNA, mt.percent)
  if(cell.annotations){
     sobj <- celltype.annotation(sobj, celltype.min, output)
  }

  return(sobj)
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
DEGage_multitest <- function(counts, group, perm.preprocess = TRUE,
                             gene.filter.threshold = 1, nperms = 2000,
                             nsubsample = NA, perm.pval = 0.1, ncores = 4,
                             maxiter = 100, mean.ratio = 1.4, subsampled.k = T,
                             writing.dir = NULL){

  #add QC's to check for number of group levels
  x <- run.multitest(counts, group,
                     perm.preprocess,
                     gene.filter.threshold,
                     nperms, nsubsample,
                     perm.pval, ncores, writing.dir)
  return(x)
}



#'@title DEGage_complete
#'@description
#'Performs cell type annotations and differential expression analysis for scRNA-seq count data
#'
#'@param input either a dataframe or matrix containing counts, or a path to an mtx directory
#'@param dir.type the type of directory that a path, if input, is directing to. Default is mtx.
#'@param celltype.min the minimum number of cells of a given cell type to count as valid.
#'@param min.nFeatureRNA The minimum number of genes that must be present in a cell
#'@param max.nFeatureRNA The maxiumum number of genes that can be present in a cell
#'@param mt.percent A value between 0 and 1 indicating the maximum proporition of mitochondrial genes allowed in a cell
#'@param cell.annotations A boolean indicating whether or not to perform celltype annotations.
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
#'@return If writing.dir is NULL, a list of dataframes containing the results for each comparison is returned. If a writing directory is provided, nothing is returned.
#'@export
DEGage_complete <- function(input, dir.type = 'mtx',
                            min.nFeatureRNA = 200, max.nFeatureRNA = 8000,
                            mt.percent = .2, cell.annotations = TRUE,
                            celltype.min = 20,  perm.preprocess = TRUE,
                            gene.filter.threshold = 1, nperms = 2000,
                            nsubsample = NA, perm.pval = 0.1, ncores = 4,
                            maxiter = 100, mean.ratio = 1.4, subsampled.k = T,
                            writing.dir = NULL){

  if(is.null(writing.dir)){
    warning("No writing directory was provided. Providing one is strongly recommended")
  }else{
    if(!dir.exists(writing.dir)){
      stop("Provided writing directory does not exist")
    }
  }

  processed <- DEGage_preprocess(input, dir.type, min.nFeatureRNA,
                                 max.nFeatureRNA, mt.percent,
                                 cell.annotations, celltype.min,
                                 output = "df")

  group <- factor(as.character(processed[[2]]))
  counts <- processed[[1]]
  rm(processed)

  print.labels(group)
  x <- run.multitest(counts, group,
                  perm.preprocess,
                  gene.filter.threshold,
                  nperms, nsubsample,
                  perm.pval, ncores,writing.dir)

  return(x)
}

