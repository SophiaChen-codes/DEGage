//Alicia Petrany, draft as of 11/13/2022
//Some functions in this file exist as facets of DR. Che's original work, and
//therefore are not used in the main DEGage function in R
//The following functions are not called in DEGage: hypergeometric, DEGage_var,
//DEGage_pdf_staticK, DEGage_cdf


#include <iostream>
#include <cmath>
#include <Rcpp.h>
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
