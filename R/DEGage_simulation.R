#Alicia Petrany, 5/18/2023

#'@title QC_Simulation
#'@description
#'Performs quality control on user input for DEGage_simulation
#'@param ngenes is the number of genes to be simulated total
#'@param ndegs is the number of genes to be differentially expressed
#'@param cellgroups is factor with integers corresponding to cell conditions
#'@param lfc is either an integer or vector with a length of ndegs which dictates log fold change value to adjust means for differentially expressed genes
#'@param prop.zeros The proportion of 0's to introduce to simulate dropouts. Must be a value between 0 and 1.
#'@param seed sets a seed for random generation. If user does not input one, a random seed is generated
#'@param ncores the number of cores to allocated for parallel computing
QC_simulation <- function(ngenes, ndegs, ngroup1, ngroup2, lfc, min.prop.zeros, max.prop.zeros, seed, ncores, dispersions, means){
  if(!inherits(ngenes, "numeric")) stop("'ngenes' must be numeric")
  if(!inherits(ndegs, "numeric")) stop("'ndegs' must be numeric")
  if(!inherits(ngroup1, "numeric")) stop("'ngroup1' must be numeric")
  if(!inherits(ngroup2, "numeric")) stop("'ngroup2' must be numeric")
  if(!inherits(min.prop.zeros, "numeric")) stop("'min.prop.zeros' must be numeric")
  if(!inherits(max.prop.zeros, "numeric")) stop("'max.prop.zeros' must be numeric")
  if(length(lfc) != 1  & length(lfc) != ndegs) stop("lfc must be an integer, or a vector the length of ndegs")
  if(min.prop.zeros > 1 & min.prop.zeros < 0) stop("min.prop.zeros must be between 0 and 1")
  if(max.prop.zeros > 1 & max.prop.zeros < 0) stop("max.prop.zeros must be between 0 and 1")
  if(length(dispersions) != 1  & length(dispersions) != ndegs & !is.null(dispersions)) stop("dispersions must be an integer, or a vector the length of ndegs")
  if(length(means) != 1  & length(means) != ndegs & !is.null(means)) stop("means must be an integer, or a vector the length of ndegs")
}

#'@title sim.core.config
#'@description sets up cores for parallelization
#'@param ncores the number of cores to use
sim.core.config <- function(ncores){
  cl.cores <- detectCores()
  if(cl.cores < ncores){
    stop("'ncores' is greater than the number of available cores. Change ncores to a smaller value")
  }
  cl <- makeCluster(getOption("cl.cores", ncores))  #allocate cores to use
  clusterExport(cl = cl,
                c("gen_rownames",
                  "gen_colnames"),
                envir = environment())
  return(cl)
}

#'@title gen_rownames
#'@description
#'A helper function that generates row names for simulation.
#' @param ngenes the number of total genes
#' @param ndegs the number of desired differentially expressed genes
#' @param ncores the number of cores to allocated for parallel computing
gen_rownames <- function(ngenes, ndegs, cl){

  Degnames <- parSapply(cl = cl,
                        X = c(1:ndegs),
                        FUN = function(x) paste("DEG",x, sep = ""))
  nondegnames <- parSapply(cl = cl,
                           X = c((ndegs+1): ngenes),
                           FUN = function(x) paste("GENE",x, sep = "" ))
  rowtitles <- append(Degnames, nondegnames)
  return(rowtitles)
}

#'@title gen_colnames
#'@description
#' A helper function that generates a name for each cell for simulation. Names indicate which condition/group each cell belongs to
#' @param ncells the number of cells to be generated
#' @param ngroup1 the number of cells in group one
#' @param ncores the number of cores to allocated for parallel computing
gen_colnames <- function(ncells, ngroup1, cl){

  group1.titles <- parSapply(cl = cl,
                             X = 1:ngroup1,
                             FUN = function(x) paste("Cell",x,".Group1", sep = "" ))
  group2.titles <- parSapply(cl = cl,
                             X = (ngroup1+1):ncells,
                             FUN = function(x) paste("Cell",x,".Group2", sep = "" ))
  celltitles <- append(group1.titles, group2.titles)
  return(celltitles)
}

#'@title get.expression.counts
#'@description
#'A helper functions that simulates counts along an NB distribution
#'@param packed.string A string of parameters separated by underscores.
#'The first value is the number of cells, the second is the mean, and
#'the third is the dipsersion
#'
get.expression.counts <- function(packed.string){
  unpack = as.numeric(strsplit(packed.string, "_")[[1]])
  x <- rnbinom(n = unpack[1], size = unpack[3]  , mu = unpack[2])
  return(x)
}

#'@title gen_colnames
#'@description
#' A helper function that introduces dropouts at a predetermined volume
#' @param finaldf simulated count matrix
#' @param ncells the number of cells simulated
#' @param ndropouts the number of dropouts to be introduced
introduce_dropouts <- function(finaldf, ncells, min.drops, max.drops){
  ndropouts = floor(ncells*runif(1, min.drops, max.drops))
  dropout.index <- sample(c(1:ncells),ndropouts)
  finaldf[dropout.index] = 0
  return(finaldf)
}

#'@title DEGage_Simulation
#'@description
#'Simulates counts. DEG's are simulated under two conditions
#'@param ngenes The number of genes to be simulated total
#'@param ndegs The number of genes to be differentially expressed
#'@param ngroup1 The number of cells in the first group
#'@param ngroup2 The number of cells in the second group
#'@param lfc Either an integer or vector with a length of ndegs which dictates log fold change value to adjust means for differentially expressed genes
#'@param min.prop.zeros The minimum proportion of dropouts that are to be introduced for a gene
#'@param max.prop.zeros The maximum proportion of dropouts that are to be introduced for a gene
#'@param seed Sets a seed for random generation. If user does not input one, a random seed is generated
#'@param ncores The number of cores to allocated for parallel computing
#'@param dispersions An optional integer or vector the length of ngenes+ndegs that sets genewise dispersion values
#'@param means An optional integer or vector the length of ngenes+ndegs that sets the means for genewise distributions                         
#'@export
DEGage_Simulation <- function(ngenes, 
                              ndegs, 
                              ngroup1, 
                              ngroup2, 
                              lfc = 2, 
                              min.prop.zeros = 0.1, 
                              max.prop.zeros = 0.5, 
                              seed = NULL, 
                              ncores = 4, 
                              dispersions = NULL, 
                              means = NULL){

  QC_simulation(ngenes, ndegs, ngroup1, ngroup2,
                lfc, min.prop.zeros, max.prop.zeros,
                seed, ncores, dispersions, means)

  #Setting the seed and other useful variables
  if(is.null(seed)){
    seed <- sample(1:100000, 1)
  }
  set.seed(seed)

  ncells = ngroup1+ngroup2

  #Detecting cores for parallel computing
  cl <- sim.core.config(ncores)

  #Name generation
  rowtitles <- gen_rownames(ngenes+ndegs, ndegs, cl)
  celltitles <- gen_colnames(ncells, ngroup1, cl)


  #Parameter Sampling. means are sampled from a gamma distrbution,
  #then multiplied by a random scale factor.
  if(is.null(dispersions)){
    size.list <- rgamma(ngenes+ndegs,1)
  }else{
    size.list <- dispersions
    if(length(size.list == 1)){
      size.list <- rep(size.list, ndegs+ngenes)
     }
  }

  if(is.null(means)){
     scale.factor.list <- sample(0:100,(ngenes+ndegs), replace = TRUE)
     mean.list <- floor(rgamma(ngenes+ndegs, 1) * scale.factor.list)
  }else{
    if(length(means) == 1){
      mean.list <- rep(means, ngenes+ndegs)
    }else{
      mean.list <- means
    }
  }


  #Group one simulation: simulates counts based off a negative binomial distrubution with predetermined means.
  #Dispersions are sampled from a gamma distrbution
  packed.means.disps <- paste(ngroup1,
                              mean.list[1:ndegs],
                              size.list[1:ndegs],
                              sep = "_")

  DEGgroup1df <- parSapply(cl = cl,
            X = packed.means.disps,
            FUN = get.expression.counts)

  DEGgroup1df <- data.frame(t(DEGgroup1df))

  #Group two simulation: simulated in the same manner as group one simulation, except means are multiplied by specified
  #LFC value
  Group2means.adj <- mean.list[1:ndegs] * (2^lfc)

  packed.means.disps <- paste(ngroup2,
                              Group2means.adj,
                              size.list[1:ndegs],
                              sep = "_")

  DEGgroup2df <- parSapply(cl = cl,
                           X = packed.means.disps,
                           FUN = get.expression.counts)

  DEGgroup2df <- data.frame(t(DEGgroup2df))

  finaldf <- cbind(DEGgroup1df, DEGgroup2df)

  #Non DEG simulation: genecounts are simulated as one large group with no change in means.

  packed.means.disps <- paste(ngroup2 + ngroup1,
                              mean.list[(ndegs + 1):(ngenes+ndegs)],
                              size.list[(ndegs + 1):(ngenes+ndegs)],
                              sep = "_")
  nondegcounts <- parSapply(cl = cl,
                           X = packed.means.disps,
                           FUN = get.expression.counts)

  nondegcounts <- data.frame(t(nondegcounts))

  #Output format
  colnames(nondegcounts) <- celltitles
  colnames(finaldf) <- celltitles
  finaldf <- rbind(finaldf, nondegcounts)
  rownames(finaldf) <- rowtitles

  #dropout management
  finaldf <- parApply(cl=cl, X = finaldf,
                      MARGIN = 1,
                      FUN = introduce_dropouts,
                      ncells = ncells,
                      min.drops = min.prop.zeros,
                      max.drops = max.prop.zeros)
  finaldf <- as.data.frame(t(as.data.frame(finaldf)))

  #Closing parallel computing cluster
  stopCluster(cl)

  return(finaldf)
}


