#Alicia Petrany, draft as of 4.20.2023
#Contains all helper function for DEGage

#'@title parse.input
#'@description
#'parses input of the user based on the data type of the input
#'
#'@param counts user input for 'counts'
#'@returns dataframe containing counts
parse.input <- function(counts){
  dtype <- class(counts)
  if(dtype != "character"){
    if(dtype == 'data.frame') return(counts)
    if(dtype == 'Seurat') return(counts@assays$RNA@counts)
    if(dtype == "SingleCellExperiment") return(counts(counts))
  } else{
    return(as.data.frame(Read10X(counts)))
  }
  stop("'counts' is not of a supported data type. See the documentation at https://rpubs.com/aliciaprowan/1041880 for supported data types")
}

#'@title QC
#'@description
#'Checks user input to ensure DEGage can run properly
#'
#'@param counts A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@return Does not return anything
QC <- function(counts, group, gene.filter.threshold, nperms){
  if(!inherits(group, "factor")) stop("'group' must be a factor. Try factor(group)")
  if(length(levels(group)) < 2) stop("'group' has less than 2 levels. DEGage requires 'group' to have 2 levels")
  if(gene.filter.threshold >1 | gene.filter.threshold < 0) {
    stop("gene.filter.threshold must be a value between 0 and 1")
  }
  if(nperms <0){
    stop("'nperms' must be greater than zero")
  }
}

#'@title subsampling
#'@description
#'Samples counts if the size of either group is greater than nsubsample
#'@param counts A dataframe containing all read counts
#'@param group A factor containing column annotations
#'@param nsubsample Number of cells to subsample
#'@return A list containing subsampled counts at index 1 and a revised group factor at index 2
subsampling <- function(counts, group, nsubsample){

  subsampling1 = FALSE
  subsampling2 = FALSE

  if(is.na(nsubsample)){
    nsubsample = 125 + ceiling(ncol(counts)*0.1)
  }

  #slices dataframes for each group individually if they have
  #more columns than nsubsample
  if(length(which(group == levels(group)[1])) > nsubsample){
    g1 <- subsampling_helper(counts, group, nsubsample, 1)
    subsampling1 = TRUE
  }

  if(length(which(group == levels(group)[2])) > nsubsample){
    g2 <- subsampling_helper(counts, group, nsubsample, 2)
    subsampling2 = TRUE
  }

  #slicing subsampled output depending on whether only one or both
  #groups required subsampling. Also fixes group
  if(subsampling1 & subsampling2){
    counts <- cbind(g1, g2)
    group <- factor(c(rep(1, nsubsample), rep(2, nsubsample)))
  }else if(subsampling1 & !subsampling2){
    counts <- cbind(g1, counts[, which(group == levels(group)[2])])
    group <- factor(c(rep(1,nsubsample),
                      rep(2, length(which(group == levels(group)[2])))))

  }else if(subsampling2 & !subsampling1){
    counts <- cbind(counts[, which(group == levels(group)[1])], g2)
    group <- factor(c( rep(1, length(which(group == levels(group)[1]))),
                       rep(2,nsubsample)))
  }

  return(list(counts, group))
}

#'@title subsampling_helper
#'@description
#'slices counts by randomly sampling columns
#'@param counts A dataframe containing all read counts
#'@param group A factor containing column annotations
#'@param nsubsample Number of cells to subsample
#'@param n group number. Either 1 or 2.
#'@return A dataframe containing containing sliced counts for one group
subsampling_helper <- function(counts, group, nsubsample, n){

  g <- counts[, which(group == levels(group)[n])]
  g <- g[,sample(1:ncol(g), nsubsample)]

  return(g)
}

#'@title core.config
#'@description
#'configures cores for parallel computing
#'@param ncores Number of cores to use
#'@param counts A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@return A cluster object to use for future parallel computing
core.config <- function(ncores, nperms, counts, group){
  cl.cores <- detectCores() #detect cores

  if(cl.cores < ncores){
    stop("'ncores' is greater than the number of available cores. Change ncores to a smaller value")
  }

  cl <- makeCluster(getOption("cl.cores", ncores));  #allocate cores to use
  clusterExport(cl = cl,
                c("counts",
                  "permtest_facilitator",
                  "group",
                  "nperms",
                  "glm.nb",
                  "zeroinfl",
                  "permtest"),
                envir = environment())
  return(cl)
}


#'@title run.permtest
#'@description
#'The first function called to initiate the permutation test, does not call c++ functions
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@param perm.pval Pvalues
#'@param nperms An integer which determines how many iterations are performed for the permutation test.
#'@param cl parallel computing cluster information
#'@return A dataframe containing gene names with their corresponding pvalues
run.permtest <- function(counts, group, perm.pval, nperms, cl){
  message("Running Permutation Test")
  permresults <- as.numeric(parApply(cl = cl,
                                     X = counts,
                                     MARGIN = 1,
                                     FUN = permtest_facilitator,
                                     group = group,
                                     nperms = nperms))
  othertail <- 1-permresults
  permresults[which(othertail < permresults)] <- othertail[which(othertail < permresults)]
  permresultdf <- data.frame("pval" = permresults, "gene" = rownames(counts))

  return(permresultdf)
}

#'@title permtest_facilitator
#'@description
#'A helper function that calls upon the permutation test written in C++, runs a single gene at a time
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@param nperms An integer which determines how many iterations are performed for the permutation test.
#'@return a pvalue fora single gene
#'@import Rcpp
permtest_facilitator <- function(counts, group, nperms){
  counts <- as.numeric(counts)
  x <- permtest(counts,
                as.numeric(group),
                levels(group),
                nperms)
  return(x)
}


#'@title run.NB.fitting
#'@description
#'calls upon NB_model_fitting in DEGage.R to perform regression on each gene
#'
#'@param counts A dataframe containing counts for each gene
#'@param group A factor that contains grouping information for the counts
#'@param cl parallel computing cluster information
#'@return dataframe with regression parameters for every gene
run.NB.fitting <- function(counts, group, cl){
  message("NB Fitting")
  counts.holder <- split(counts, seq(nrow(counts)))
  counts.holder <-lapply(X = counts.holder,
                         FUN = function(x) as.numeric(x))
  outputdf <- parLapplyLB(cl=cl,
                          X = counts.holder,
                          fun = NB_model_fitting,
                          group = group)
  outputdf <- as.data.frame(matrix(unlist(outputdf),
                                   nrow = length(outputdf),
                                   byrow = TRUE))

  rownames(outputdf) <- rownames(counts)
  colnames(outputdf) <- c("r1", "p1","mu1", "r2", "p2", "mu2", "base.mean")
  return(outputdf)
}


#'@title run.NB.fitting
#'@description
#'calls upon NB_model_fitting in DEGage.R to perform regression on each gene
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@param outputdf A dataframe containing regression parameters and pvalues
#'@param cl parallel computing cluster information
#'@return a dataframe containing all regression paramters and pvalues
run.ZINB.refitting <- function(counts,group, outputdf, cl){
  message("ZINB Fitting")
  temp.holder <- split(counts, seq(nrow(counts)))
  temp.holder <-  lapply(X = temp.holder, FUN = function(x) as.numeric(x))
  temp.outputdf <- parLapplyLB( cl=cl,
                                X = temp.holder,
                                fun =ZINB_model_fitting,
                                group = group)
  temp.outputdf <- as.data.frame(matrix(unlist(temp.outputdf),
                                        nrow = length(temp.outputdf),
                                        byrow = TRUE))
  colnames(temp.outputdf) <- c("z.r1", "z.p1", "z.r2", "z.p2")

  return(cbind(outputdf, temp.outputdf))
}

#'@title fix_infs
#'@description
#'Fixes genes that could not be fit with DOTNB because of limitations with regression
#'or hypergeometric function calcualtions.
#'@param outputdf A dataframe containing regression parameters and pvalues
#'@param mean.ratio A positive value indicating the minimum mean ratio for an gene with an incalculable pvalue to be considered differentially expressed
fix_infs <- function(outputdf, mean.ratio){
  outputdf[ (is.infinite(outputdf$pval) | is.na(outputdf$pval)) &
              (outputdf$mu1 != 0 & outputdf$mu2 != 0 ) &
              (outputdf$mu1/outputdf$mu2 > mean.ratio | outputdf$mu2/outputdf$mu1 > mean.ratio), ]$FDR <-
    rep(0, nrow(outputdf[ (is.infinite(outputdf$pval) | is.na(outputdf$pval)) &
                            (outputdf$mu1 != 0 & outputdf$mu2 != 0 ) &
                            (outputdf$mu1/outputdf$mu2 > mean.ratio | outputdf$mu2/outputdf$mu1 > mean.ratio), ] ))

  outputdf[ (is.infinite(outputdf$pval) | is.na(outputdf$pval)) &
              (outputdf$mu1 != 0 & outputdf$mu2 != 0 ) &
              (outputdf$mu1/outputdf$mu2 > mean.ratio | outputdf$mu2/outputdf$mu1 > mean.ratio), ]$pval <-
    rep(0, nrow(outputdf[ (is.infinite(outputdf$pval) | is.na(outputdf$pval)) &
                            (outputdf$mu1 != 0 & outputdf$mu2 != 0 ) &
                            (outputdf$mu1/outputdf$mu2 > mean.ratio | outputdf$mu2/outputdf$mu1 > mean.ratio), ] ))
  return(outputdf)
}


#'@title organize_output
#'@description
#'formats output for user
#'@param outputdf df with all information to be output to user
#'@return An outputdf with the columns in a different order
organize_output <- function(outputdf){
  pval <- outputdf$pval
  permPvals <- outputdf$permPvals
  FDR <- outputdf$FDR

  outputdf <- outputdf[,1:11]

  # outputdf$fit.method <- fit.method
  outputdf$permPvals <- permPvals
  outputdf$pval <- pval
  outputdf$FDR <- FDR

  return(outputdf)
}



