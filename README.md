# DEGage
This package allows for the differential expression analysis of two groups of count data based on a difference of two negative binomial distributions. 
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
