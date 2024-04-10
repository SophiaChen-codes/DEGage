#Alicia Petrany, draft as of 4.5.2024
#This file contains mathematical components of the model

#'@title NB_model_fitting
#'@description
#'A helper function that performs genewise standard negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import stats
#'@import MASS
#'@return A dataframe with a single row containing NB regression parameters for a gene
NB_model_fitting <- function(counts, group){
  #Isolate counts as vector
  groupone <- as.numeric(counts[which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[which(group == levels(group)[2])])

  if(sum(groupone)==0){
    #add a single, very small psuedocount so dist can still exist
    groupone[1] = groupone[1] + 0.00000001
  }
  if(sum(grouptwo)==0){
    #add a single, very small psuedocount so dist can still exist
    grouptwo[1] = grouptwo[1] + 0.00000001
  }

  df1 <- data.frame(counts = as.numeric(groupone))
  df2 <- data.frame(counts = as.numeric(grouptwo))

  #tryCatch must be used during regression, because in cases where the data is not negative binomial,
  #MASS will throw an error. NA's are returned a regression parameters in this case.
  r1 <- tryCatch({r1<-summary(glm.nb(counts ~1, data = df1))$theta},
                 error = function(e){return(NA)})
  r2 <- tryCatch({r2<-summary(glm.nb(counts ~1, data = df2))$theta},
                 error = function(e){return(NA)})

  mu1 <- mean(groupone)
  mu2 <- mean(grouptwo)

  if(is.na(r1)){
    r1 = NA
    p1 = NA
  }else{
    p1 <- r1/(r1+mu1)
  }

  if(is.na(r2)){
    r2 = NA
    p2 = NA
  }else{
    p2 <- r2/(r2+mu2)
  }

  output <- data.frame(r1 = r1,
                       p1 = p1,
                       mu1 = mu1,
                       r2 = r2,
                       p2 = p2,
                       mu2 = mu2,
                       basemean = mean(c(groupone, grouptwo)),
                       row.names = rownames(counts))
  return(output)
}

#'@title ZINB_model_fitting
#'@description
#'A helper function that performs genewise zero-inflated negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import pcsl
#'@import stats
#'@return A dataframe with a single row containing regression parameters for a gene
ZINB_model_fitting <- function(counts, group){
  #Identical to NB_model_fitting, except PSCL's zeroinfl function is used for regression
  groupone <- as.numeric(counts[which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[which(group == levels(group)[2])])

  if(sum(groupone)==0){
    #add a single, very small psuedocount so dist can still exist
    groupone[1] = groupone[1] + 0.00000001
  }
  if(sum(grouptwo)==0){
    #add a single, very small psuedocount so dist can still exist
    grouptwo[1] = grouptwo[1] + 0.00000001
  }

  #if no zeros in row, introduce a random dropout so zinb can be fit
  if(!(0 %in% groupone)){
    groupone[sample(c(1:length(groupone)), 1)] <- 0
  }
  if(!(0 %in% grouptwo)){
    grouptwo[sample(c(1:length(grouptwo)), 1)] <- 0
  }

  df1 <- data.frame(counts = as.numeric(ceiling(groupone)))
  df2 <- data.frame(counts = as.numeric(ceiling(grouptwo)))

  r1 = tryCatch({r1<-zeroinfl(counts ~1, data = df1, dist = "negbin")$theta},
                error = function(e){return(NA)})
  r2 = tryCatch({r2<-zeroinfl(counts ~1, data = df2, dist = "negbin")$theta},
                error = function(e){return(NA)})

  mu1 = mean(groupone)
  mu2 = mean(grouptwo)

  if(is.na(r1)){
    r1 = NA
    p1 = NA
  }else{
    p1 <- r1/(r1+mu1)
  }

  if(is.na(r2)){
    r2 = NA
    p2 = NA
  }else{
    p2 <- r2/(r2+mu2)
  }

  output <- data.frame(r1 = r1,
                       p1 = p1,
                       r2 = r2,
                       p2 = p2,
                       row.names = rownames(counts))
  return(output)
}


#'@title calculate_k_random
#'@description
#'Calculate k based on the random assignment protocol
#'
#'@param row counts for a single gene
#'@param group A factor that contains grouping information for the counts
#'@return k for each gene
calculate_k_random <- function(row, group){
  c2 <- row[group == levels(group)[1]]
  c1 <- row[group == levels(group)[2]]
  swap = FALSE
  if(length(c2) > length(c1)){
    temp <- c1
    c1 <- c2
    c2 <- temp
    swap = TRUE
  }
  if(length(c1) != length(c2)){
    c2 <- c(sample(c2, length(c2)),
            sample(c2, length(c1) - length(c2), replace = TRUE))
  }
  #shuffle
     c1 <- sample(c1, length(c1))
     c2 <- sample(c2, length(c2))
  # calculate k
  if(swap){
    return(mean(c2-c1))
  }else{
    return(mean(c1-c2))
 }
}

#'@title calculate_k_random
#'@description
#'Calculate k based on the subsampling protocol
#'
#'@param row counts for a single gene
#'@param group A factor that contains grouping information for the counts
#'@return k for each gene
calculate_k_subsample <- function(row, group){
  meanvec <- c()
    c2 <- row[group == levels(group)[1]]
    c1 <- row[group == levels(group)[2]]
    swap = FALSE
    if(length(c2) > length(c1)){
      temp <- c1
      c1 <- c2
      c2 <- temp
      swap = TRUE
    }
    if(length(c1) != length(c2)){
      c1 <- sample(c1, length(c2))
    }
    #shuffle
    c1 <- sample(c1, length(c1))
    c2 <- sample(c2, length(c2))
    # calculate k
    if(swap){
      meanvec <- c(meanvec, mean(c2-c1))
    }else{
     meanvec <- c(meanvec, mean(c1-c2))
   }
  return(mean(meanvec))
}

