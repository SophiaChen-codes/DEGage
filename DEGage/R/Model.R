#Alicia Petrany, draft as of 4.20.2023
#This file contains mathematical components of the model

#' NB_model_fitting
#'
#'A helper function that performs genewise standard negative binomial regression
#'
#'@param counts A sliced dataframe that contains the counts for an individual gene.
#'@param group A factor that contains grouping information for the counts
#'@import stats
#'@import MASS
#'@return A dataframe with a single row containing NB regression parameters for a gene
NB_model_fitting <- function(counts, group){

  #Isolate counts as vector
  print(rownames(counts))
  groupone <- as.numeric(counts[which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[which(group == levels(group)[2])])

  if((sum(groupone)==0)| (sum(grouptwo) == 0)){
    output <- data.frame(r1 = NA,
                         p1 = NA,
                         mu1 = NA,
                         r2 = NA,
                         p2 = NA,
                         mu2 = NA,
                         row.names = rownames(counts))
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

  if(is.na(r1) | is.na(r2)){
    return(data.frame(r1 = NA,
                      p1 = NA,
                      mu1 = NA,
                      r2 = NA,
                      p2 = NA,
                      mu2 = NA,
                      row.names = rownames(counts)))
  }

  mu1 <- mean(groupone)
  p1 <- r1/(r1+mu1)

  mu2 <- mean(grouptwo)
  p2 <- r2/(r2+mu2)

  output <- data.frame(r1 = r1,
                       p1 = p1,
                       mu1 = mu1,
                       r2 = r2,
                       p2 = p2,
                       mu2 = mu2,
                       row.names = rownames(counts))
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
#'@return A dataframe with a single row containing regression parameters for a gene
ZINB_model_fitting <- function(counts, group){

  #Identical to NB_model_fitting, except PSCL's zeroinfl function is used for regression
  groupone <- as.numeric(counts[which(group == levels(group)[1])])
  grouptwo <- as.numeric(counts[which(group == levels(group)[2])])

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

  if(is.na(r1) | is.na(r2)){
    return(data.frame(r1 = NA,
                      p1 = NA,
                      mu1 = NA,
                      r2 = NA,
                      p2 = NA,
                      mu2 = NA,
                      row.names = rownames(counts)))
  }

  mu1 <- mean(groupone)
  p1 <- r1/(r1+mu1)

  mu2 <- mean(grouptwo)
  p2 <- r2/(r2+mu2)

  output <- data.frame(r1 = r1,
                       p1 = p1,
                       mu1 = mu1,
                       r2 = r2,
                       p2 = p2,
                       mu2 = mu2,
                       row.names = rownames(counts))
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
#'@return A pvalue for a single gene
#'@import hypergeo

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
    pval = (p1^r1)*(p2^r2)*((gamma(r1+k)/gamma(r1))/gamma(1+k))*(q1^k)*genhypergeo(c(r1+k,r2),k+1,q1*q2,
                                                                                   maxiter = 100000,
                                                                                   polynomial = TRUE)
  }else{
    nk=-k
    pval = (p1^r1)*(p2^r2)*((gamma(r2+nk)/gamma(r2))/gamma(1+nk))*(q2^nk)*genhypergeo(c(r2+nk,r1),nk+1,q1*q2,
                                                                                      maxiter = 100000,
                                                                                      polynomial = TRUE)
  }
  return(pval)
}

#' DEGage_mean
#'
#'Calculates mean of distrubtion
#'
#'@param r1 r parameter for the first NB distribution
#'@param p1 p parameter for the first NB distribution
#'@param r2 r parameter for the second NB distrbution
#'@param p2 p parameter for the second NB distribution
#'@return a mean for a single gene
DEGage_mean <- function(r1, p1, r2, p2) {
  q1 = 1 - p1
  q2 = 1 - p2
  dmean = ((r1 * q1) / p1) - (r2 * q2 / p2)
  return(dmean)
}

