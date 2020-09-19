# IMIX
R package for integrative genomics analysis using a multivariate mixture model framework.

**Reference**

* Ziqiao Wang and Peng Wei. IMIX: A multivariate mixture model approach to integrative analysis of multiple types of omics data. bioRxiv 2020.06.23.167312; doi: https://doi.org/10.1101/2020.06.23.167312

### Installation
Install the IMIX package via github
```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ziqiaow/IMIX",build_vignettes=TRUE) 
```

### To get started with the package
Load the package and open the package vignette:
```
library("IMIX")
vignette("IMIX")
```
There will be examples for all the functions in IMIX package in [vignette](vignettes/IMIX.html). If you have any question, use help().

Have fun with IMIX :smiley:

