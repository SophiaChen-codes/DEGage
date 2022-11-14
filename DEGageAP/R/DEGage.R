#Alicia Petrany, draft as of 11/13/2022
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
#'@param na.rm a logical indicating whether or not NAs should be removed
#'@import Rcpp
#'@export
DEGage <- function(counts, group, perm.preprocess = TRUE, gene.filter.threshold = 1, nperms = 2000){

  if(!require("MASS")){
    suppressPackageStartupMessages(install.packages("MASS"))
    library("MASS")
  }
  if(!require("pscl")){
    install.packages("pscl")
    suppressPackageStartupMessages(library("pscl"))
  }
  if(!require("hypergeo")){
    install.packages("hypergeo")
    suppressPackageStartupMessages(library("hypergeo"))
  }
  counts <- counts[rowSums(counts == 0)/ncol(counts) < gene.filter.threshold,]
  outputdf <- data.frame()

  #Running gene permutation test
  permresults <- apply(counts, MARGIN = 1, FUN = permtest_facilitator, group = group, nperms = nperms, simplify = TRUE)
  permresultdf <- data.frame("pval" = as.numeric(permresults), "gene" = rownames(counts))
  if(perm.preprocess){
    counts <- counts[permresultdf$pval  < .1,]
  }

  #Performs genewise negative binomial regression, generates the df that is to be output with regression parameters
  outputdf <- apply(counts, MARGIN = 1, FUN = trad_model_fitting, group = group)
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
    fillerdf <- zeroinfl_model_fitting(counts[i,], group)
    if(is.na(fillerdf$r1)){
      fillerdf$pval <- NA
    }
    fillerdf$pval <- pdf_facilitator(fillerdf)
    outputdf[i,] <- as.numeric(fillerdf)
  }

  #If perm.preprocess is false, replaces NA's with permutation test pvalue,
  #or the maxiumum between the DEGage and the Permutation test pvalues are selected
  if(!perm.preprocess){
    for(i in 1:nrow(outputdf)){
      if( (is.na(outputdf[i,]$pval)) | (outputdf[i,]$pval < permresultdf[i,]$pval)){
        #outputdf[i,]$pval = permresultdf[i,]$pval
      }
    }
  }

  #P value correction for multiple tests
  outputdf$FDR <- p.adjust(outputdf$pval, method = "fdr")

  return(outputdf)
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

#' trad_model_fitting
#'
#'A helper function that performs genewise standard negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import stats
#'@import MASS
#'
trad_model_fitting <- function(counts, group){
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

#'zeroinfl_model_fitting
#'
#'A helper function that performs genewise zero-inflated negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import pcsl
#'@import stats
#'
zeroinfl_model_fitting <- function(counts, group){
  #Identical to trad_model_fitting, except PSCL's zeroinfl function is used for regression
  groupone <- as.numeric(counts[,which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[,which(group == levels(group)[2])])

  if((sum(groupone)==0)| (sum(grouptwo) == 0)| !(0 %in% groupone) | !(0 %in% grouptwo)){
    output <- data.frame(r1 = NA, p1 = NA, mu1 = NA, r2 = NA, p2 = NA, mu2 = NA, row.names = rownames(counts))
    return(output)
  }

  df1 <- data.frame(counts = as.numeric(ceiling(groupone)))
  df2 <- data.frame(counts = as.numeric(ceiling(grouptwo)))

  r1 = tryCatch({r1<-zeroinfl(counts ~1, data = df1, dist = "negbin")$theta},
                error = function(e){return(NA)})
  r2 = tryCatch({r2<-zeroinfl(counts ~1, data = df2, dist = "negbin")$theta},
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


#'gen_rownames
#'
#'A helper function that generates row names for simulation.
#' @param ngenes the number of total genes
#' @param ndegs the number of desired differentially expressed genes
gen_rownames <- function(ngenes, ndegs){
  Degnames <- sapply(X = c(1:ndegs), FUN = function(x) paste("DEG",x, sep = ""))
  nondegnames <- sapply(X = c((ndegs+1): ngenes), FUN = function(x) paste("GENE",x, sep = "" ))
  rowtitles <- append(Degnames, nondegnames)
  return(rowtitles)
}


#' gene.colnames
#'
#' A helper function that generates a name for each cell for simulation. Names indicate which condition/group each cell belongs to
#' @param ncells the number of cells to be generated
#' @param ngroup1 the number of cells in group one
gen_colnames <- function(ncells, ngroup1){
  group1.titles <- sapply(X = 1:ngroup1, FUN = function(x)paste("Cell",x,".Group1", sep = "" ))
  group2.titles <- sapply(X = (ngroup1+1):ncells, FUN = function(x)paste("Cell",x,".Group2", sep = "" ))
  celltitles <- append(group1.titles, group2.titles)
  return(celltitles)
}


#'DEGage_Simulation
#'
#'Simulates counts. DEG's are simulated under two conditions
#'@param ngenes is the number of genes to be simulated total
#'@param ndegs is the number of genes to be differentially expressed
#'@param cellgroups is factor with integers corresponding to cell conditions
#'@param lfc is either an integer or vector with a length of ndegs which dictates log fold change value to adjust means for differentially expressed genes
#'@param seed sets a seed for random generation. If user does not input one, a random seed is generated.
#'@import stats
#'@export

DEGage_Simulation <- function(ngenes, ndegs, cellgroups, lfc = 1, prop.zeros = .3, seed = NULL){

   if(!require("stats")){
    install.pacakges("stats")
    suppressPackageStartupMessages(library("stats"))
   }

  #Setting the seed and other useful variables
  if(is.null(seed)){
    seed <- sample(1:100000, 1)
  }
  set.seed(seed)

  ncells = length(cellgroups)
  ngroup1 <- sum(cellgroups == levels(cellgroups)[1])
  ngroup2 <- sum(cellgroups == levels(cellgroups)[2])
  finaldf <- data.frame()

  #Name generation
  rowtitles <- gen_rownames(ngenes+ndegs, ndegs)
  celltitles <- gen_colnames(ncells, ngroup1)

  #Parameter Sampling. means are sampled from a gamma distrbution,
  #then multiplied by a random scale factor.
  scale.factor.list <- sample(0:100,(ngenes+ndegs), replace = TRUE)
  mean.list <- floor(rgamma(ngenes+ndegs, 1) * scale.factor.list)
  size.list <- rgamma(ngenes+ndegs,1)

  #DEG simulation.
  DEGmeans <- mean.list[1:ndegs]

  #Group one simulation: simulates counts based off a negative binomial distrubution with predetermined means.
  #Dispersions are sampled from a gamma distrbution
  DEGgroup1list <- lapply(DEGmeans, FUN = function(x) rnbinom(n = ngroup1, size = size.list, mu = x))
  DEGgroup1df <- as.data.frame(matrix(unlist(DEGgroup1list), nrow = length(DEGgroup1list), byrow = TRUE))

  #Group two simulation: simulated in the same manner as group one simulation, except means are multiplied by specified
  #LFC value
  Group2means.adj <- DEGmeans * (2^lfc)
  DEGgroup2list <- lapply(Group2means.adj, FUN = function(x) rnbinom(n = ngroup2, size = rgamma(1,1), mu = x))
  DEGgroup2df <- as.data.frame(matrix(unlist(DEGgroup2list), nrow = length(DEGgroup2list), byrow = TRUE))

  finaldf <- cbind(DEGgroup1df, DEGgroup2df)


  #Non DEG simulation: genecounts are simulated as one large group with no change in means.
  nondeglist <- lapply(mean.list[(ndegs+1):length(mean.list)], FUN = function(x) rnbinom(n = ncells, size =  runif(1, min = .01, max = 2), mu = x))
  nondegdf <- as.data.frame(matrix(unlist(nondeglist), nrow = length(nondeglist), byrow = TRUE))

  #Output format
  colnames(nondegdf) <- celltitles
  colnames(finaldf) <- celltitles
  finaldf <- rbind(finaldf, nondegdf)
  rownames(finaldf) <- rowtitles

  #dropout management
  ndropouts = ncells*prop.zeros
  for(i in 1:ndegs){
    dropout.index <- sample(c(1:ncells),ndropouts)
    finaldf[i, dropout.index] = 0
  }

  return(finaldf)
}


