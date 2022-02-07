# Dynamic spectrum of ectopic lymphoid B cell activation in the RA synovium characterized by NR4A nuclear receptor expression

Preprint is [here](https://www.biorxiv.org/content/10.1101/2021.05.14.443150v2)



## Organization and use

### [synovium_sc](synovium_sc/)
code and primary data for expression and BCR for single cells from peripheral blood and synovium from RA134, RA172, RA195, RA221, e.g. figure 1A-1D, 2A-2D.
Uses R 3.5.1

### [reanalysis](reanalysis/) 
Code to infer cell types on external data sets e.g. figure 1E.
Uses R 3.6.1

### [amp_bulk](amp_bulk/)
differential expression of NR4A in low-input bulk data, e.g. figure 1F. 
Uses R 4.1.


## Data

[Raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196150) from the 10X run in this manuscript is on GEO.  Most the batteries are included to run these analysis here.

Exceptions: AMP data for the amp_bulk and reanalysis workflows will need to be downloaded from IMMPort [1](https://www.immport.org/shared/study/SDY998), [2](https://www.immport.org/shared/study/SDY1299) and/or GEO, since we did not have permission to directly (re)-share.

## R packages

To install necessary packages to reproduce an analysis, change to a subdirectories, and run:
```r
install.packages('renv')
renv::restore()
```

You will need to be running the same version of R as used when the renv lockfiles were written...

