# DEGage
This package allows for the differential expression analysis of two groups of scRNA-seq count data. It employs a novel family of discrete distributions for describing the difference of two NB distributions (named DOTNB). DEGage take the raw counts of scRNA-seq as inputs, and thus avoid introducing artificially bias in normalization steps in current methods. A workflow is shown as follows.

![DEGage Workflow](/Fig1_DEGage_Workflow.png)

## Install
To install DEGage, run the following:
```
library(devtools)
install_github("chenyongrowan/DEGage")
```
## Documentation
For more detailed documentation, see:    
https://rpubs.com/aliciaprowan/1043456

DEGage contains 5 main functions: DEGage(), DEGage_multitest(), DEGage_preprocess(), DEGage_complete(), and DEGage_Simulation().   

### DEGage
DEGage performs pairwise differential analysis on NGS count data. The input is typically a dataframe where columns contain samples and rows contain genes. Counts do not need to be normalized prior to use with DEGage. 
```
DEGage(counts, group, perm.preprocess = TRUE,
       gene.filter.threshold = 1, nperms = 2000,
       nsubsample = NA, perm.pval = 0.1, ncores = 4,
       maxiter = 1000, mean.ratio = 1.4, subsampled.k = T)
```

### DEGage_multitest
Similar to DEGage, except it can perform an indefinite number of pairwise comparisons. 
```
DEGage_multitest(counts, group, perm.preprocess = FALSE,
                 gene.filter.threshold = 1, nperms = 2000,
                 nsubsample = NA, perm.pval = 0.1, ncores = 4,
                 writing.dir = NULL)
```

### DEGage_preprocess
Performs a simple scRNA-seq count pre-processing workflow with Seurat and generates automatic cell type annotations with SingleR. It is primarily used in DEGage_complete; however, it is available for use individually.  
```
DEGage_preprocess(input, dir.type = 'mtx', min.nFeatureRNA = 200,
                  max.nFeatureRNA = 8000, mt.percent = .2,
                  cell.annotations = TRUE,
                  celltype.min = 20, output = "Seurat")
```

### DEGage_complete
Takes unprocessed scRNA-seq counts as an input, automatically processes and generates cell type annotations for them, then performs differential expression analysis on the groups identified through cell-type clustering.    
```
DEGage_complete(input, dir.type = 'mtx',
                min.nFeatureRNA = 200, max.nFeatureRNA = 8000, 
                mt.percent = .2, cell.annotations = TRUE,
                celltype.min = 20, perm.preprocess = FALSE,
                gene.filter.threshold = 1, nperms = 2000, 
                nsubsample = NA, perm.pval = 0.1,
                ncores = 4, writing.dir = NULL)
```

### DEGage_Simulation
Generates simplified simulated scRNA-Seq counts following an NB distribution with pre-defined proportions of dropouts.   
```
DEGage_Simulation(ngenes, ndegs, cellgroups, lfc = 1, 
                  prop.zeros = .3, seed = NULL, ncores = 4)
```

## Example Usage
In this section, we will detail how to use DEGage functions

First, we will simulate a small data frame of counts to pass through DEGage() using DEGage simulate. It has 100 cells total with 50 in each condition, as well as 1000 genes with 100 degs: 
```
library(DEGage)
cellgroups <- factor(c(rep(1,50), rep(2,50)))
df <- DEGage_Simulation(ngenes = 1000, ndegs = 100, cellgroups = cellgroups)
```

Next, we will pass these counts through DEGage:
```
results <- DEGage(counts = df, group = cellgroups)
```

To test DEGage_multitest, we will simulate a second dataframe of counts, merge them together, and pass them through DEGage_multitest(): 
```
cellgroups2 <- factor( c( rep(3,50), rep(4,50) ) )
df2 <- DEGage_Simulation(ngenes = 1000, ndegs = 100, cellgroups = cellgroups2)

df <- cbind(df, df2)
cellgroups <- factor(c(cellgroups, cellgroups2))

multitest.results <- DEGage_multitest(df, cellgroups)
```
## Citation
Please cite the following article if you use DEGage in your research:

Petrany A., Zhang S. and Chen, Y. A General Probabilistic Model for Comparative Analysis of Next-generation Sequencing data. Under review. 
