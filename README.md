# DEGage
This package allows for the differential expression analysis of two groups of scRNA-seq count data. It employs a novel family of discrete distributions for describing the difference of two NB distributions (named DOTNB). DEGage take the raw counts of scRNA-seq as inputs, and thus avoid introducing artificially bias in normalization steps in current methods. A workflow is shown as follows.

![DEGage Workflow](/DEGage/DEGage_Workflow.png)

## Install
To install DEGage, copy and paste the following into your R terminal:
```
library(devtools)
install_github("chenyongrowan/DEGage")
```
## Documentation
For detailed documentation, see:    
https://rpubs.com/aliciaprowan/1043456

## Citation
Please cite the following article if you use DEGage in your research:

Petrany A., Zhang S. and Chen, Y. DEGage: a General Model-based Method for Detecting Differentially Expressed Genes from scRNA-seq Data. To be submitted. 
