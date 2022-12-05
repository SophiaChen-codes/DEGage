//Alicia Petrany, draft as of 11/13/2022
//Some functions in this file exist as facets of DR. Che's original work, and
//therefore are not used in the main DEGage function in R
//The following functions are not called in DEGage: hypergeometric, DEGage_var,
//DEGage_pdf_staticK, DEGage_cdf


#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <omp.h>. ## this is used for parallel computing
using namespace Rcpp;

//'permtest
//'
//'Performs genewise permutation tests to calculate pvalues for prefiltering.
//'
//'@param genecount An R Numeric vector containing counts for a single gene
//'@param group An R Numeric vector containing group info, typically two integers that correspond to groups
//'@param levels The levels aspect of an R factor for the above group parameter
//'@param nperms The number of permutations to perform
// [[Rcpp::export]]
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
    // following statement used for parallel computing.
    #pragma omp parallel for
    for (int i = 0; i < nperms; i++) {
        intermediate = sample(genecount, genecount.length(), false);
        g1 = intermediate[Range(0, (genecount.length() / 2) - 1)];
        g2 = intermediate[Range(genecount.length() / 2, genecount.length() - 1)];
        perm = abs((mean(g2) - mean(g1)));
        if (perm > trueval) {
            larger++;
        }
    }

    return(larger / nperms);
}

double hypergeometric(double a, double b, double c, double x)
{
    const double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    int n = 1;

    while (abs(term) > TOLERANCE)
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / n;
        value += term;
    }
    return value;
}

// [[Rcpp::export]]
double DEGage_mean(double r1, double p1, double r2, double p2) {
    double q1 = 1 - p1, q2 = 1 - p2;
    double dmean = ((r1 * q1) / p1) - (r2 * q2 / p2);
    return dmean;
}
// [[Rcpp::export]]
double DEGage_var(double r1, double p1, double r2, double p2) {
    double q1 = 1 - p1, q2 = 1 - p2;
    double dvar = (r1 * q1) / pow(p1,2) + (r2 * q2) / pow(p2,2);
    return dvar;
}

double DEGage_pdf_staticK(double r1, double p1, double r2, double p2, int k) {
    double q1 = 1 - p1, q2 = 1 - p2;
    double prob;

    if (k > 0) {
        prob = pow(p1, r1) * pow(p2, r2) * ((tgamma(r1 + k) / tgamma(r1)) / tgamma(1 + k)) * pow(q1, k) * hypergeometric(r1 + k, r2, k + 1, q1 * q2);
    }
    else {
        k = -k;
        prob = pow(p1, r1) * pow(p2, r2) * ((tgamma(r2 + k) / tgamma(r2)) / tgamma(1 + k)) * pow(q2, k) * hypergeometric(r2 + k, r1, k + 1, q1 * q2);
    }
    return prob;
}

// [[Rcpp::export]]
double DEGage_cdf(double r1, double p1, double r2, double p2, int k) {
    float stop_pvalue = 1.0 * pow(10, -20);
    double pdfval = DEGage_pdf_staticK(r1, p1, r2, p2, k);
    double cvalue = 0;

    if (k < DEGage_mean(r1, p1, r2, p2)){
        if (pdfval < stop_pvalue) {
            cvalue = cvalue + pdfval;
            k = k-1;
        }
        else {
            while (pdfval > stop_pvalue) {
                pdfval = DEGage_pdf_staticK(r1, p1, r2, p2, k);
                cvalue = cvalue + pdfval;
                k = k - 1;
            }
        }
    }
    else {
        if (pdfval < stop_pvalue) {
            cvalue = 1;
        }
        else {
            while (DEGage_pdf_staticK(r1, p1, r2, p2, k) > stop_pvalue) {
                pdfval = DEGage_pdf_staticK(r1, p1, r2, p2, k);
                cvalue = cvalue + pdfval;
                k = k - 1;
            }
        }
    }
    return(cvalue);
}
