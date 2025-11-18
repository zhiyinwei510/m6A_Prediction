# m6APrediction

## Purpose

`m6APrediction` is an R package for predicting RNA m6A modification
sites using a Random Forest model.  
It uses sequence and structure features such as GC content, RNA type,
RNA region, exon length, distance to junction, evolutionary
conservation, and DNA 5mer to estimate the probability and
classification of m6A sites.

## ROC Curve
![ROC Curve](https://github.com/zhiyinwei510/m6APrediction/blob/main/inst/images/roc_curve.png)

## PRC Curve
![PRC Curve](https://github.com/zhiyinwei510/m6APrediction/blob/main/inst/images/prc_curve.png)

------------------------------------------------------------------------

## Installation

Install the package directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("zhiyinwei510/m6APrediction")
```

## Example Usage

### Run multiple m6A prediction examples

``` r
rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
prediction_multiple(rf_fit, example_data, 0.6)
```

### Run single m6A prediction example

``` r
rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
prediction_single(rf_fit, 0.6, "mRNA", "CDS", 12, 5, 0.8, "ATCGAT", 0.5)
```
