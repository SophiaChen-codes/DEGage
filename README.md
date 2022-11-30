# DEGage
This package allows for the differential expression analysis of two groups of scRNA-seq count data. It employs a novel family of discrete distributions for describing the difference of two NB distributions (named DOTNB). DEGage take the raw counts of scRNA-seq as inputs, and thus avoid introducing artificially bias in normalization steps in current methods. A workflow is shown as follows.

![alt text](https://github.com/chenyongrowan/DEGage/tree/main/DEGageAP/DEGage_Workflow.tif)

 
## Requirements
```
pscl 1.5.5
MASS 7.3-58.1
stats 4.2.2
hypergeo 1.2-13
```
DEGage will install these dependencies automatically upon the first time running it if they are not already installed
## Main Functions
#### DEGage
```
DEGage(counts, group, perm.preprocess = TRUE, gene.filter.threshold = 1, nperms = 2000)
```
Tests for differential expression between two groups of count data  
##### Parameters:  
- `counts`: A data frame where samples correspond to columns and genes correspond to rows   
- `group`: A factor of a numeric vector that indicates which samples belong to which pairwise condition. An example is shown below under the Example Usage section  
- `perm.preprocess`: A boolean indicating whether a permuation test is used to prefilter genes  
- `gene.filter.threshold`: A value between 0-1 indicating the maximum proportion of zero counts a gene can have before being filtered out. A value of 1 means only genes with all zero counts are filtered out, where as a value of 0.5 would filter out genes where half the counts are zeros  
- `nperms`: The number of permutations performed during the permutation test  
#### DEGage_Simulation
```
DEgage_Simulation(ngenes, ndegs, cellgroups, lfc = 1, prop.zeros = .3, seed = NULL)
```
Simulates counts following a negative binomial distribution with predefined levels of differential expression and zero count proportions
##### Parameters
- `ngenes`: The number of equivalently expressed genes to simulate  
- `ndegs`: The number of differentially expressed genes to simulate  
- `cellgroups`: A factor structured similarly to the group parameter from DEGage. It is a factor of a numeric vector that contains integers corresponding to the conditions of each cell. This vector also inherently defines the number of cells to generate  
- `lfc`: The log fold change that differentially expressed genes are up/down regulated by. Can be an integer for uniform expression changes across all degs, or a vector the same length as ndegs that specifies lcf values for each differentially expressed gene  
- `prop.zeros`: A value between 0-1 that specifies the proportion of dropouts to be added to each gene  
- `seed`: Specify a seed for random generation to generate a reproducible matrix  
## Example Usage
### Differential Expression Analysis with DEGage
As an example case, we will use a scRNA-seq dataset from Rao-Ruiz et. al that contains 38 neurons from the denate gyrus of *Mus musculus*. For this study, we would like to examine differential expression between dVenus+/- cells from the fear condtioned (FC) mice.  
First, we will load in the DEGage library, as well as stringi for some basic string manipulation
```
library(DEGageAP)
library(stringi)
```
Next, we will retrieve the Rao-Ruiz dataset from the GEO
```
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129024/suppl/GSE129024_counts_per_gene_sample.txt.gz"
tmp <- tempfile()
download.file(url, tmp)
counts <- read.delim(gzfile(tmp), sep = "\t")
```
Next, we will isolate the desired conditions and compile them into a new dataframe
```
#Isolating FC Neurons
FCcounts <- counts[,substr(colnames(counts), 1, 2)=="FC"]
#dVenus+ Neurons
dvenus.pos <- FCcounts[,substr(colnames(FCcounts), 4,4) == "G"]
#dVenus- Neurons
dvenus.neg <- FCcounts[,substr(colnames(FCcounts),4,4)=="N"]
FCcounts <- cbind(dvenus.neg, dvenus.pos)
```
Now, we will generate the factor for DEGage's `group` parameter
```
groups <- c(rep(1, ncol(dvenus.neg)),rep(2,ncol(dvenus.pos)))
groups <- factor(groups)
```
Next we will run DEGage
```
DEGage.Results <- DEGage(counts = FCcounts, group = groups)
```
DEGage.Results now contains a dataframe containng the NB regression outputs (mu, r, and p), a p-value, and an FDR for each gene.  
##### Some Important Notes:
- While `perm.preprocess` = TRUE, the dataframe returned from DEGage will likely contain signficantly fewer genes than the input dataframe. This is because the permutation test eliminates these genes prior to regression. This step increases DEGage's specficity. However, if you would like a p-value for each gene, you should set `perm.preprocess` to FALSE. 
### Count simulation with DEGage
The `cellgroups` parameter is defined in the same manner as the `group` parameter for DEGage. If you needed to produce a count data frame with 100 cells equally divided into 2 conditions, 18000 equivalently expressed genes, and 2000 differentially expressed genes, it would be done as follows: 
```
library(DEGageAP)

Cellgroups <- factor(c(rep(1,50),rep(2,50)))
counts <- DEGage_Simulation(ngenes = 18000, ndegs = 2000, cellgroups = Cellgroups)
```
