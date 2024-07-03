# DEGage
This package allows for the differential expression analysis of two groups of scRNA-seq count data. It employs a novel family of discrete distributions for describing the difference of two NB distributions (named DOTNB). DEGage take the raw counts of scRNA-seq as inputs, and thus avoid introducing artificially bias in normalization steps in current methods. A workflow is shown as follows.

![DEGage Workflow](/Fig1_DEGage_Workflow.png)

## Installation
To install DEGage, run the following:
```
library(devtools)
install_github("chenyongrowan/DEGage")
```
## Dependencies
DEgage requires the following dependencies:
```
MASS >= 7.3-59
pscl >= 1.5.5
hypergeo >= 1.2-13
stats >= 4.2.2
doParallel >= 1.0.17
parallel >= 4.3.0
foreach >=1.5.2
```
A recent MASS update made it incompatible with R versions below 4.4. To install a legacy version of MASS, go to https://cran.r-project.org/src/contrib/Archive/MASS/ and download an archived version that is compatible with your current R installation. Then, run the following in an R terminal: 
```
install.pacakges("path/to/MASS_SOURCE.tar.gz", repos=NULL, type="source")
```
## Differential Expression Functions
For more detailed documentation and tutorials, please refer to RPubs by RStudio:    
https://rpubs.com/aliciaprowan/1043456

### DEGage
DEGage performs pairwise differential analysis on NGS count data. The input is typically a dataframe where columns contain samples and rows contain genes. DEGage equips with two sampling strategies, random assignment and subsampling, for handling imbalanced scRNA-seq datasets.
```
DEGage(counts,
     group,
     perm.preprocess = TRUE,
     gene.filter.threshold = 0.9,
     nperms = 2000,
     nsubsample = NA,
     perm.pval = 0.1,
     ncores = 4,
     maxiter = 100,
     mean.ratio = 4,
     subsampled.k = T)
```

### DEGage_multitest
Similar to DEGage, except it can perform an indefinite number of pairwise comparisons. 
```
DEGage_multitest(counts,
               group,
               perm.preprocess = TRUE,
               gene.filter.threshold = 0.9,
               nperms = 2000,
               nsubsample = NA,
               perm.pval = 0.1,
               ncores = 4,
               maxiter = 100,
               mean.ratio = 4,
               subsampled.k = T,
               writing.dir = NULL)
```

### DEGage_preprocess
An optional function that performs basic prefiltering on a dataframe of read counts prior to calling DEGage()
```
DEGage_preprocess(counts,
                min.transcripts.per.gene = 200,
                max.transcripts.per.gene = 10000,
                min.transcripts.per.cell = 500,
                max.transcripts.per.cell = 20000,
                min.cells.per.gene = 3,
                max.dropout.prop = 0.8,
                max.mt.percent = 0.2,
                mt.prefix = 'MT-',
                max.rRna.percent = 0.2,
                rRna.prefix = 'rr')
```


### DEGage_Simulation
Generates simplified simulated scRNA-Seq counts following an NB distribution with pre-defined proportions of dropouts.   
```
DEGage_Simulation(ngenes,
       ndegs,
       ngroup1,
       ngroup2,
       lfc = 1,
       min.prop.zeros = 0.1,
       max.prop.zeros = 0.5,
       seed = NULL,
       ncores = 4,
       dispersions = NULL,
       means = NULL)
```

## Example Usage
In this section, we will detail how to use DEGage functions

First, we will simulate a small data frame of counts to pass through DEGage() using DEGage simulate. It has 100 cells total with 50 in each condition, as well as 1000 genes with 100 degs: 
```
library(DEGage)
df <- DEGage_Simulation(ngenes = 1000, ndegs = 100, ngroup1 = 50, ngroup2 = 50)
```

Next, we will pass these counts through DEGage:
```
cellgroups <- factor(c(rep(1,50), rep(2,50)))
results <- DEGage(counts = df, group = cellgroups)
```

To test DEGage_multitest, we will simulate a second dataframe of counts, merge them together, and pass them through DEGage_multitest(): 
```
df2 <- DEGage_Simulation(ngenes = 1000, ndegs = 100, ngroup1 = 50, ngroup2 = 50)

df <- cbind(df, df2)
cellgroups2 <- factor(c(rep(3,50), rep(4,50)))
cellgroups <- factor(c(cellgroups, cellgroups2))

multitest.results <- DEGage_multitest(df, cellgroups)
```
## Citation
Please cite the following article if you use DEGage in your research:

Petrany A., Chen R., Zhang S. and Chen, Y. "Theoretical Framework for the Difference of Two Negative Binomial Distributions and its Application in Comparative Analysis of Sequencing Data". Genome Research, Under Review.
