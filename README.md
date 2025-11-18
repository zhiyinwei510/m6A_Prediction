# m6APrediction

## Purpose

`m6APrediction` is an R package for predicting N6-methyladenosine (m6A) RNA modification sites using machine learning models such as Random Forest. It provides functions to encode DNA sequences, perform multiple-sample prediction, and make single-sample predictions with interpretable probability outputs.

## ROC Curve
![ROC Curve](https://github.com/zhiyinwei510/m6A_Prediction/blob/main/inst/images/roc_curve.png)

## PRC Curve
![PRC Curve](https://github.com/zhiyinwei510/m6A_Prediction/blob/main/inst/images/prc_curve.png)

------------------------------------------------------------------------

## Installation

Install the package directly from GitHub:

``` r
#install.packages("remotes")
remotes::install_github("zhiyinwei510/m6A_Prediction")
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
