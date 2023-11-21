#Alicia Petrany, draft as of 4.20.2023
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

  if((sum(groupone)==0)| (sum(grouptwo) == 0)){
    output <- data.frame(r1 = NA,
                         p1 = NA,
                         mu1 = mean(groupone),
                         r2 = NA,
                         p2 = NA,
                         mu2 = mean(grouptwo),
                         base.mean = mean(c(groupone, grouptwo)),
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

  if((sum(groupone)==0)| (sum(grouptwo) == 0)| !(0 %in% groupone) | !(0 %in% grouptwo)){
    output <- data.frame(r1 = NA,
                         p1 = NA,
                         r2 = NA,
                         p2 = NA,
                         row.names = rownames(counts))
    return(output)
  }

  df1 <- data.frame(counts = as.numeric(ceiling(groupone)))
  df2 <- data.frame(counts = as.numeric(ceiling(grouptwo)))

  r1 = tryCatch({r1<-zeroinfl(counts ~1, data = df1, dist = "negbin")$theta},
                error = function(e){return(NA)})
  r2 = tryCatch({r2<-zeroinfl(counts ~1, data = df2, dist = "negbin")$theta},
                error = function(e){return(NA)})

  mu1 = mean(groupone)
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
                       r2 = r2,
                       p2 = p2,
                       row.names = rownames(counts))
  return(output)
}

