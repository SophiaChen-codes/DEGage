#Alicia Petrany, draft as of 11/29/2022
#Install on machine with devtools::install()
#'DEGage
#'
#'Tests for pairwise differential expression between two count distributions
#'
#'@param input A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param perm.preprocess A logical indicating whether or not to perform the permutation prefiltering step.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
#'@import Rcpp
#'@export
DEGage <- function(counts, group, perm.preprocess = TRUE, gene.filter.threshold = 1, nperms = 2000){

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
  ##Quality control of user input
  QC(counts, group, gene.filter.threshold ,nperms)

  counts <- counts[rowSums(counts == 0)/ncol(counts) < gene.filter.threshold,]
  outputdf <- data.frame()

  #Running gene permutation test
  if(perm.preprocess){
    permresults <- as.numeric(apply(counts, MARGIN = 1, FUN = permtest_facilitator, group = group, nperms = nperms, simplify = TRUE))
    othertail <- 1-permresults
    permresults[which(othertail < permresults)] <- othertail[which(othertail < permresults)]
    permresultdf <- data.frame("pval" = permresults, "gene" = rownames(counts))
    counts <- counts[permresultdf$pval  < 0.05,]
  }

  #Performs genewise negative binomial regression, generates the df that is to be output with regression parameters
  outputdf <- apply(counts, MARGIN = 1, FUN = NB_model_fitting, group = group)
  outputdf <- as.data.frame(matrix(unlist(outputdf), nrow = length(outputdf), byrow = TRUE))
  rownames(outputdf) <- rownames(counts)
  colnames(outputdf) <- c("r1", "p1","mu1", "r2", "p2", "mu2")
  message("Finished Fitting")

  #Calculates pvals based off of regression parameters calculated above, adds new column to outputdf containing pvals
  outputdf$pval <- pdf_facilitator(outputdf)

  message("refitting")

  #Refits genes with NA's after the first regression with zero inflated negative binomial regression
  NA.indices <- which(is.na(outputdf$pval))
  for(i in NA.indices){
    fillerdf <- ZINB_model_fitting(counts[i,], group)
    if(is.na(fillerdf$r1)){
      fillerdf$pval <- NA
    }
    fillerdf$pval <- pdf_facilitator(fillerdf)
    outputdf[i,] <- as.numeric(fillerdf)
  }

  #P value correction for multiple tests
  outputdf$FDR <- p.adjust(outputdf$pval, method = "fdr")

  return(outputdf)
}

#'QC
#'
#'Checks user input to ensure DEGage can run properly
#'
#'@param input A dataframe with cells as columns and genes as rows.
#'@param group A factor which assigns conditions to cells as in Deseq2.
#'@param gene.filter.threshold A value between 0-1 which represents the maximum proportion of zeros a gene can have before being filtered out.
#'@param nperms An integer greater than 0 that indicates how many permutations will be carried out during the permutation test
QC <- function(counts, group, gene.filter.threshold, nperms){
  if(ncol(counts)!= length(group)) stop("The length of 'group' and the number of columns in counts must be equal")
  if(class(group) != "factor") stop("'group' must be a factor. See ... for example")
  if(length(levels(group)) > 2) stop("'group' has more than 2 levels. DEGage only tests for pairwise condtions")
  if(length(levels(group)) < 2) stop("'group' has less than 2 levels. DEGage requires 'group' to have 2 levels")
  if(class(counts) != "data.frame") stop("'counts' must be a data.frame")
  if(gene.filter.threshold >1 | gene.filter.threshold < 0) stop("gene.filter.threshold must be a value between 0 and 1")
  if(nperms <0) stop("'nperms' must be greater than zero")
}

#' permtest_facilitator
#'
#'A helper function that calls upon the permutation test written in C++
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@param nperms An integer which determines how many iterations are performed for the permutation test.
#'
#'@import Rcpp
#'
permtest_facilitator <- function(counts, group, nperms){
  counts <- as.numeric(counts)
  x <- permtest(counts, as.numeric(group), levels(group), nperms)
  return(x)
}

#' NB_model_fitting
#'
#'A helper function that performs genewise standard negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import stats
#'@import MASS
#'
NB_model_fitting <- function(counts, group){
  #Isolate counts as vector
  groupone <- as.numeric(counts[which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[which(group == levels(group)[2])])

  if((sum(groupone)==0)| (sum(grouptwo) == 0)){
    output <- data.frame(r1 = NA, p1 = NA, mu1 = NA, r2 = NA, p2 = NA, mu2 = NA, row.names = rownames(counts))
    return(output)
  }

  df1 <- data.frame(counts = as.numeric(groupone))
  df2 <- data.frame(counts = as.numeric(grouptwo))

  #tryCatch must be used during regression, because in cases where the data is not negative binomial,
  #MASS will throw an error. NA's are returned a regression parameters in this case.
  r1 <- tryCatch({r1<-summary(glm.nb(counts ~1, data = df1))$theta},
                 error = function(e){return(NA)})

  r2 <- tryCatch({r2<-summary(glm.nb(counts ~1, data = df2))$theta},
                 error = function(e){return(NA)})

  if(is.na(r1) || is.na(r2)){
    return(data.frame(r1 = NA, p1 = NA, mu1 = NA, r2 = NA, p2 = NA, mu2 = NA, row.names = rownames(counts)))
  }

  mu1 <- mean(groupone)
  p1 <- r1/(r1+mu1)

  mu2 <- mean(grouptwo)
  p2 <- r2/(r2+mu2)

  output <- data.frame(r1 = r1, p1 = p1, mu1 = mu1, r2 = r2, p2 = p2, mu2 = mu2, row.names = rownames(counts))
  return(output)
}

#'ZINB_model_fitting
#'
#'A helper function that performs genewise zero-inflated negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import pcsl
#'@import stats
#'
ZINB_model_fitting <- function(counts, group){
  #Identical to NB_model_fitting, except PSCL's ZINB function is used for regression
  groupone <- as.numeric(counts[,which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[,which(group == levels(group)[2])])

  if((sum(groupone)==0)| (sum(grouptwo) == 0)| !(0 %in% groupone) | !(0 %in% grouptwo)){
    output <- data.frame(r1 = NA, p1 = NA, mu1 = NA, r2 = NA, p2 = NA, mu2 = NA, row.names = rownames(counts))
    return(output)
  }

  df1 <- data.frame(counts = as.numeric(ceiling(groupone)))
  df2 <- data.frame(counts = as.numeric(ceiling(grouptwo)))

  r1 = tryCatch({r1<-ZINB(counts ~1, data = df1, dist = "negbin")$theta},
                error = function(e){return(NA)})
  r2 = tryCatch({r2<-ZINB(counts ~1, data = df2, dist = "negbin")$theta},
                error = function(e){return(NA)})
  if(is.na(r1) || is.na(r2)){
    return(data.frame(r1 = NA, p1 = NA, mu1 = NA, r2 = NA, p2 = NA, mu2 = NA, row.names = rownames(counts)))
  }


  mu1 <- mean(groupone)
  p1 <- r1/(r1+mu1)


  mu2 <- mean(grouptwo)
  p2 <- r2/(r2+mu2)

  output <- data.frame(r1 = r1, p1 = p1, mu1 = mu1, r2 = r2, p2 = p2, mu2 = mu2, row.names = rownames(counts))
  return(output)
}

#' DEGage_pdf
#'
#'Calculates pvalues from regression parameters
#'
#'@param r1 r parameter for the first NB distribution
#'@param p1 p parameter for the first NB distribution
#'@param r2 r parameter for the second NB distrbution
#'@param p2 p parameter for the second NB distribution
#'@param k Observed number of differences
#'
#'@import hypergeo
#'
DEGage_pdf<-function(r1,p1,r2,p2,k){

  q1 = 1-p1
  q2 = 1-p2
  #Calculations with |k| >10 provide a level of precision that is not necessary,
  #and High k's throw errors in the genhypergeo function
  if(k >10){
    k = 10
  }else if(k < -10){
    k = -10
  }
  if(k>0){
    pval = (p1^r1)*(p2^r2)*((gamma(r1+k)/gamma(r1))/gamma(1+k))*(q1^k)*genhypergeo(c(r1+k,r2),k+1,q1*q2, maxiter = 100000, polynomial = TRUE)
  }else{
    nk=-k
    pval = (p1^r1)*(p2^r2)*((gamma(r2+nk)/gamma(r2))/gamma(1+nk))*(q2^nk)*genhypergeo(c(r2+nk,r1),nk+1,q1*q2, maxiter = 100000, polynomial = TRUE)
  }
  return(pval)
}

#'pdf_facilitator
#'
#'A helper function that calls upon DEGage_pdf
#'
#'@param df The data frame containing regression parameters for every gene
#'@import Rcpp
#'
pdf_facilitator <- function(df){
  pdf.vals <- vector(mode = "numeric")
  for(i in 1:nrow(df)){
    if(is.na(df$r1[i])){
      pdf.vals <- c(pdf.vals, NA)
    }else{
      r1 <- df$r1[i]
      p1 <- df$p1[i]
      r2 <- df$r2[i]
      p2 <- df$p2[i]
      k <- DEGage_mean(r1, p1, r2, p2)
      pdf.vals<- c(pdf.vals, DEGage_pdf(r1, p1, r2, p2, k))
    }
  }
  return(pdf.vals)
}

