// [[Rcpp::interfaces(r, cpp)]]

#include <math.h>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
float DEGage_pdf(double r1, double p1, double r2, double p2, double dn );

//'@name permtest
//'@title permtest
 //'@description
 //'Performs genewise permutation tests to calculate pvalues for prefiltering.
 //'
 //'@param genecount An R Numeric vector containing counts for a single gene
 //'@param group An R Numeric vector containing group info, typically two integers that correspond to groups
 //'@param levels The levels aspect of an R factor for the above group parameter
 //'@param nperms The number of permutations to perform
 //'@export
 //[[Rcpp::export]]
 double permtest(NumericVector genecount, NumericVector group, CharacterVector levels, int nperms) {

   NumericVector group1;
   NumericVector group2;

   int group1index = 0;
   int group2index = 0;

   for (int i = 0; i < genecount.length(); i++) {
     if (group[i] == atoi(levels[0])) {
       group1.insert(group1index, genecount[i]);
       group1index++;
     }
     else {
       group2.insert(group2index, genecount[i]);
       group2index++;
     }
   };

   double trueval = mean(group2) - mean(group1);

   NumericVector perms;
   NumericVector intermediate;
   NumericVector g1;
   NumericVector g2;
   double perm;
   double larger = 0.0;

   for (int i = 0; i < nperms; i++) {
     intermediate = sample(genecount, genecount.length(), false);
     g1 = intermediate[Range(0, (genecount.length() / 2) - 1)];
     g2 = intermediate[Range(genecount.length() / 2, genecount.length() - 1)];
     perm = abs((mean(g2) - mean(g1)));
     if (perm > trueval) {
       larger++;
     }
   }

  if((larger / nperms) < (1-(larger/nperms))){
    return(larger / nperms);
  }else{
    return(1-(larger/nperms));
  }
 }

 //'@name DEGage_mean
 //'@title Calculates mean of DOTNB distribution
 //'@description
 //'Calculates the mean for the difference of two negative binomial distributions
 //'@param r1 r parameter for the first group
 //'@param p1 p parameter for the first group
 //'@param r2 r parameter for the second group
 //'@param p2 p parameter for the second group
 //'@export
 //[[Rcpp::export]]
 double DEGage_mean(double r1, double p1, double r2, double p2){
   double q1 = 1 - p1;
   double q2 = 1 - p2;
   double dmean = (r1 * q1 / p1) - (r2 * q2 / p2);
   return dmean;
 }

//'@name hypergeo
//'@title hypergeo
//'@description
//'calculates an approximation for the hypergeometric function for use in DEGage_pdf. Equivalent
//'to a 2F1 hypergeometric function calcuation with series expansion.
//'@param a Upper bound 1
//'@param b Upper bound 2
//'@param c lower bound
//'@param x primary complex argument
 // Note convergence restrictions: abs(x) < 1 and c not a negative integer or zero
 double hypergeo( double a, double b, double c, double x )
 {
   int maxiter = 10000000;
   const double TOLERANCE = 1.0e-10;
   double term = a * b * x / c;
   double value = 1.0 + term;
   int n = 1;
   int iter = 0;

   while ( abs( term ) > TOLERANCE )
   {
     a++, b++, c++, n++;
     term *= a * b * x / c / n;
     value += term;
     iter++;
     if(iter > maxiter){
       break;
     }
   }
   return value;
 }

 //'@name adjust_r_lower
 //'@title adjust_r_lower
 //'@description
 //'If one r value is too large for a pvalue to be calculated, the r value is adjusted until a p-value is calculable.
 //'This function operates for upper tailed tests.
 //'@param r1 r parameter for the first group
 //'@param p1 p parameter for the first group
 //'@param r2 r parameter for the second group
 //'@param p2 p parameter for the second group
 //'@param dn The observed number of difference between the two groups
float adjust_r_upper(double r1, double p1, double r2, double p2, double dn){
  double q1 = 1-p1, q2 = 1-p2;
  float prob = 1.1;
  double increment;

  if(r1 > 1000 && r2 > 1000){
    while(prob > 1){
      increment = r1/10;
      r1 = r1 - increment;
      r2 = r2 - increment;
      if( r1 < 1 || r2 < 1){
        break;
      }
      prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r1+dn)/tgamma(r1))/tgamma(1+dn))*pow(q1,dn)*hypergeo(r1+dn, r2, dn+1, q1*q2);
    }
  }
  if(r1 > r2){
    while(prob > 1){
      increment = r1/10;
      r1 = r1 - increment;
      if( r1 < 1){
        break;
      }
      prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r1+dn)/tgamma(r1))/tgamma(1+dn))*pow(q1,dn)*hypergeo(r1+dn, r2, dn+1, q1*q2);

    }
  }else{
    while(prob > 1){
      increment = r1/10;
      r1 = r1 - increment;
      if( r1 < 1){
        break;
      }
      prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r1+dn)/tgamma(r1))/tgamma(1+dn))*pow(q1,dn)*hypergeo(r1+dn, r2, dn+1, q1*q2);

    }
  }


  if (prob > 1){
    return(NA_REAL);
  }
  return prob;
}

//'@name adjust_r_lower
//'@title adjust_r_lower
//'@description
//'If one r value is too large for a pvalue to be calculated, the r value is adjusted until a p-value is calculable.
//'Generally only necessary in cases of extreme differential expression, hypergeometric function will not converge otherwise
//'This function operates for lower tailed tests.
//'@param r1 r parameter for the first group
//'@param p1 p parameter for the first group
//'@param r2 r parameter for the second group
//'@param p2 p parameter for the second group
//'@param dn The observed number of difference between the two groups
 float adjust_r_lower(double r1, double p1, double r2, double p2, double dn){
   double q1 = 1-p1, q2 = 1-p2;
   float prob = 1.1;
   double increment;

    if(r1 > 1000 && r2 > 1000){
      while(prob > 1){
      increment = r1/10;
      r1 = r1 - increment;
      r2 = r2 - increment;
      if( r1 < 1 || r2 < 1){
        break;
      }
      prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r2+dn)/tgamma(r2))/tgamma(1+dn))*pow(q2,dn)*hypergeo(r2+dn, r1, dn+1, q1*q2);
    }
    }
   if(r1 > r2){
    while(prob > 1){
      increment = r1/10;
        r1 = r1 - increment;
       if( r1 < 1){
         break;
       }
       prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r2+dn)/tgamma(r2))/tgamma(1+dn))*pow(q2,dn)*hypergeo(r2+dn, r1, dn+1, q1*q2);

   }
   }else{
     while(prob > 1){
       increment = r1/10;
       r1 = r1 - increment;
       if( r1 < 1){
         break;
       }
       prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r2+dn)/tgamma(r2))/tgamma(1+dn))*pow(q2,dn)*hypergeo(r2+dn, r1, dn+1, q1*q2);

     }
   }

   if (prob > 1){
     return(NA_REAL);
   }
   return prob;
 }

 //'@name DEGage_pdf
 //'@title Calculates PDF of DOTNB distribution
 //'@description
 //'Calculates the pdf for the difference of two negative binomial distributions
 //' //'@param r1 r parameter for the first group
 //'@param p1 p parameter for the first group
 //'@param r2 r parameter for the second group
 //'@param p2 p parameter for the second group
 //'@param dn The observed number of difference between the two groups
 //'@export
 // [[Rcpp::export]]
float DEGage_pdf(double r1, double p1, double r2, double p2, double dn )
 { // calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
   // dn is the observed number of difference
   double q1 = 1-p1, q2 = 1-p2;
   float prob;
   if(dn>0){
     prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r1+dn)/tgamma(r1))/tgamma(1+dn))*pow(q1,dn)*hypergeo(r1+dn, r2, dn+1, q1*q2);
     if(prob > 1){
        prob = adjust_r_upper(r1, p1, r2, p2, dn);
     }
   }else{
     dn=-dn;
     prob = pow(p1,r1)*pow(p2,r2)*((tgamma(r2+dn)/tgamma(r2))/tgamma(1+dn))*pow(q2,dn)*hypergeo(r2+dn, r1, dn+1, q1*q2);
     if(prob > 1){
       prob = adjust_r_lower(r1, p1, r2, p2, dn);
     }

   }
   return prob;
 }


 //'@name DEGage_cdf
 //'@title Calculates CDF of DOTNB distribution
 //'@description
 //'Calculations the CDF for the difference of two negative binomial distributions
 //'@param r1 r parameter for the first group
 //'@param p1 p parameter for the first group
 //'@param r2 r parameter for the second group
 //'@param p2 p parameter for the second group
 //'@param k The observed number of difference between the two groups
 //'@param maxiter The maximum number of iterations to perform during cdf calculation
 //'@export
 // [[Rcpp::export]]
float DEGage_cdf(double r1, double p1, double r2, double p2, double k, int maxiter){
   float cvalue = 0;
   float pdfv = DEGage_pdf(r1, p1, r2, p2, k);
   int iter = 0;
   k = round(k);

   if(NumericVector::is_na(pdfv) || pdfv == 0){
     if(k == 0){
       return NA_REAL;
     }
     pdfv = DEGage_pdf(r1, p1, r2, p2, k);
   }

     if(k < DEGage_mean(r1, p1, r2, p2)){
       while(!NumericVector::is_na(pdfv)){
         cvalue=cvalue+pdfv;
         k = k - 1;
         pdfv=DEGage_pdf(r1,p1,r2,p2,k);
         iter = iter + 1;
         if(iter > maxiter){
           break;
         }
       }

     }else{
       while(!NumericVector::is_na(pdfv)){
         cvalue=cvalue+pdfv;
         k = k + 1;
         pdfv=DEGage_pdf(r1,p1,r2,p2,k);
         iter = iter + 1;
         if(iter > maxiter){
           break;
         }
       }
     }

   if(cvalue > 1){
      cvalue = NA_REAL;
   }

   return cvalue;
 }

 //'@name get_min
 //'@title get_min
 //'@Description
 //'Gets the minimum between pvalues calculated with NB and ZINB parameters,
 //'excludes 0 values from assessment
 //'@param pvals A NumericVector that contains NB and ZINB fitted pvalues
 double get_min(NumericVector pvals){
   double minimum = min(na_omit(pvals));
   if(minimum == 0){
     while(minimum == 0){
       pvals = pvals[!is_nan(pvals)];
       pvals = pvals[!is_infinite(pvals)];
       if(pvals.length() < 2){
         return 0;
       } else {
         pvals.erase(which_min(pvals));
         minimum = min(na_omit(pvals));
       }
     }
   }
   return minimum;
 }


 //'@name cdf_facilitator
 //'@title cdf_facilitator
 //'@description
 //'Calculates pvalues with the CDF of the DOTNB distribution
 //'
 //'@param df A Dataframe containing regression parameters
 //'@param maxiter The maxiumum number of iterations to perform for cdf
 // [[Rcpp::export]]
 NumericVector cdf_facilitator(DataFrame& df,  int maxiter){

   NumericVector r1 = df["r1"];
   NumericVector r2 = df["r2"];
   NumericVector p1 = df["p1"];
   NumericVector p2 = df["p2"];
   NumericVector mu1 = df["mu1"];
   NumericVector mu2 = df["mu2"];
   NumericVector basemean = df["base.mean"];
   NumericVector k = df["k"];

   NumericVector z_r1 = df["z.r1"];
   NumericVector z_r2 = df["z.r2"];
   NumericVector z_p1 = df["z.p1"];
   NumericVector z_p2 = df["z.p2"];

   NumericVector cdfvals;

   for(int i = 0; i < df.nrows(); i++){
         if((NumericVector::is_na(r1[i]) || NumericVector::is_na(r2[i])) &&
        (NumericVector::is_na(z_r1[i]) || NumericVector::is_na(z_r2[i])) ){
       cdfvals.push_back(NA_REAL);

     } else if(NumericVector::is_na(r1[i]) || NumericVector::is_na(r2[i])){
       cdfvals.push_back(DEGage_cdf(z_r1[i], z_p1[i], z_r2[i], z_p2[i], k[i], maxiter));
     } else if(NumericVector::is_na(z_r1[i]) || NumericVector::is_na(z_r2[i])){
       cdfvals.push_back(DEGage_cdf(r1[i], p1[i], r2[i], p2[i], k[i], maxiter));
     }else{
       NumericVector pvals = {
         DEGage_cdf(r1[i], p1[i], r2[i], p2[i], k[i], maxiter),
         DEGage_cdf(z_r1[i], z_p1[i], z_r2[i], z_p2[i], k[i], maxiter)
       };
       cdfvals.push_back(get_min(pvals));
     }

   }

   return cdfvals;
 }

