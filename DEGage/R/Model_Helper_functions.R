#Alicia Petrany, draft as of 4.20.2023
#Contains all helper function for DEGage


#'.init.dependencies
#'
#'Loads in dependenies that are not already loaded into the environment.
#'Installs them automatically if they are unavilable
#'
#'@return Does not return anything
init.dependencies <- function(){
  if(!suppressPackageStartupMessages(require("MASS"))){
    suppressPackageStartupMessages(install.packages("MASS"))
    library("MASS")
  }
  if(!suppressPackageStartupMessages(require("pscl"))){
    install.packages("pscl")
    suppressPackageStartupMessages(library("pscl"))
  }
  if(!suppressPackageStartupMessages(require("hypergeo"))){
    install.packages("hypergeo")
    suppressPackageStartupMessages(library("hypergeo"))
  }
  if(!suppressPackageStartupMessages(require("parallel"))){
    install.packages("parallel")
    suppressPackageStartupMessages(library("parallel"))
  }
  if(!suppressPackageStartupMessages(require("doParallel"))){
    install.packages("doParallel")
    suppressPackageStartupMessages(library("doParallel"))
  }
  if(!suppressPackageStartupMessages(require("stringi"))){
    install.packages("stringi")
    suppressPackageStartupMessages(library("stringi"))
  }
}

#'QC
#'
#'Checks user input to ensure DEGage can run properly
#'
#'@param input A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@return Does not return anything
QC <- function(counts, group, gene.filter.threshold, nperms){
  if(class(group) != "factor") {
    warning("'group' is not a factor. Converting to factor.")
    group <- factor(group)
  }
  if(length(levels(group)) < 2) stop("'group' has less than 2 levels. DEGage requires 'group' to have 2 levels")
  if(class(counts) != "data.frame") stop("'counts' must be a data.frame")
  if(gene.filter.threshold >1 | gene.filter.threshold < 0) {
    stop("gene.filter.threshold must be a value between 0 and 1")
  }
  if(nperms <0){
    stop("'nperms' must be greater than zero")
  }
}

#'subsampling
#'
#'Samples counts if the size of either group is greater than nsubsample
#'@param counts A dataframe containing all read counts
#'@param group A factor containing column annotations
#'@param nsubsample Number of cells to subsample
#'@return A list containing subsampled counts at index 1 and a revised group factor at index 2
subsampling <- function(counts, group, nsubsample){

  subsampling1 = FALSE
  subsampling2 = FALSE

  if(is.na(nsubsample)){
    nsubsample = 125 + ceiling(ncol(counts)*0.01)
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

#'subsampling_helper
#'
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

#'core.config
#'
#'configures cores for parallel computing
#'@param ncores Number of cores to use
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
                  "genhypergeo",
                  "permtest"),
                envir = environment())
  return(cl)
}

#'run.permtest
#'
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

#' permtest_facilitator
#'
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

#'pdf_facilitator
#'
#'A helper function that calls upon DEGage_pdf
#'
#'@param df The data frame containing regression parameters for every gene
#'@param ncores number of cores to used
#'@param cl Parallel computing cluster info
#'@import Rcpp
pdf_facilitator <- function(df, ncores, cl){
  pdf.vals <- vector( mode = "numeric")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  pdf.vals <- foreach(i = c(1:nrow(df)),
              .export = c("genhypergeo", "DEGage_mean", "DEGage_pdf"),
              .packages = c("hypergeo")) %dopar% {

                 if(is.na(df$r1[i])){
                    pdf.vals <- NA
                 }else{
                    r1 <- df$r1[i]
                    p1 <- df$p1[i]
                    r2 <- df$r2[i]
                    p2 <- df$p2[i]
                    k <- DEGage_mean(r1, p1, r2, p2)
                    pdf.vals <- DEGage_pdf(r1, p1, r2, p2, k)
                 }
  }
  stopCluster(cl)
  registerDoSEQ()
  pvals <- unlist(pdf.vals)
  return(pvals)
}

#'run.NB.fitting
#'
#'calls upon NB_model_fitting in DEGage.R to perform regression on each gene
#'
#'@param counts A dataframe containing counts for each gene
#'@param group A factor that contains grouping information for the counts
#'@param cl parallel computing cluster information
#'@return dataframe with regression parameters for every gene
run.NB.fitting <- function(counts, group, cl){
  message("Fitting Counts")
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
  colnames(outputdf) <- c("r1", "p1","mu1", "r2", "p2", "mu2")
  return(outputdf)
}

#'run.NB.fitting
#'
#'calls upon NB_model_fitting in DEGage.R to perform regression on each gene
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@param outputdf A dataframe containing regression parameters and pvalues
#'@param cl parallel computing cluster information
#'@param ncores number of cores to use in parallel computing
#'@return a dataframe containing all regression paramters and pvalues
run.ZINB.refitting <- function(counts,group, outputdf, cl, ncores){
  message("Refitting Counts")
  NA.indices <- which(is.na(outputdf$pval))
  temp <- counts[NA.indices,]
  temp.holder <- split(temp, seq(nrow(temp)))
  temp.holder <-  lapply(X = temp.holder, FUN = function(x) as.numeric(x))
  temp.outputdf <- parLapplyLB( cl=cl,
                                X = temp.holder,
                                fun =ZINB_model_fitting,
                                group = group)
  temp.outputdf <- as.data.frame(matrix(unlist(temp.outputdf),
                                        nrow = length(temp.outputdf),
                                        byrow = TRUE))
  colnames(temp.outputdf) <- c("r1", "p1","mu1", "r2", "p2", "mu2")
  temp.outputdf$pval <- pdf_facilitator(temp.outputdf,
                                        ncores,
                                        cl)
  outputdf[NA.indices,] <- temp.outputdf
  outputdf$fit.method <- rep("NB", nrow(outputdf))
  outputdf$fit.method[NA.indices] <- rep("ZINB", length(NA.indices))
  return(outputdf)
}

#'organize_output
#'
#'@param outputdf df with all information to be output to user
#'@return An outputdf with the columns in a different order
organize_output <- function(outputdf){
  pval <- outputdf$pval
  fit.method <- outputdf$fit.method
  permPvals <- outputdf$permPvals
  FDR <- outputdf$FDR

  outputdf <- outputdf[,1:6]

  outputdf$fit.method <- fit.method
  outputdf$permPvals <- permPvals
  outputdf$pval <- pval
  outputdf$FDR <- FDR

  return(outputdf)
}
