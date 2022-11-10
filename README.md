# IMIX
R package for integrative genomics analysis using a multivariate mixture model framework. This package is also available in R/Cran: https://cran.rstudio.com/web/packages/IMIX/index.html.

**Reference**

* Ziqiao Wang and Peng Wei. IMIX: A multivariate mixture model approach to integrative analysis of multiple types of omics data. Bioinformatics. 2021 Apr 1;36(22-23):5439-5447. doi: 10.1093/bioinformatics/btaa1001. PMID: 33258948; PMCID: PMC8016490.

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
There will be examples for all the functions in IMIX package in [vignette](vignettes/IMIX.pdf). If you have any question, use help().

Have fun with IMIX :smiley:

