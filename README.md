# slmrf

This repository contains the code used for the simulation study in the manuscript "Comparing Spatial Regression to Random Forests for Large Environmental Data Sets."  The repository is set up as an R package that contains functions that can simulate spatially autocorrelated data, estimate spatial regression models, and make kriging predictions.  The simulation code (R Markdown file) and data are located in the `\inst` folder of the package directory.

To install:

* using the `devtools` package: `install_github("ericwfox/slmrf")`

In addition to the simulation code, the cleaned MMI response data and StreamCat covariates, which were used to fit the spatial regression and random forest models in the paper, can be accessed in this R package.  To load the data sets run the following commands:

`data("mmi_data")`

`data("streamcat_preds")`

For users that are not familiar with R, the data are also available in CSV format in the `\inst` folder.  



