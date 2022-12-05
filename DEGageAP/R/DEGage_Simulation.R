#'gen_rownames
#'
#'A helper function that generates row names for simulation.
#' @param ngenes the number of total genes
#' @param ndegs the number of desired differentially expressed genes

library(parallel)
cl.cores <- detectCores()   #detect available cores
cl <- makeCluster(getOption("cl.cores", 4));  #allocate cores to use
### change the following sapply, lapply to parSapply and parLapply accordingly.

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
