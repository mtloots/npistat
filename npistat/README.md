
# npistat

[![R-CMD-check](https://github.com/mtloots/npistat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/npistat/actions/workflows/R-CMD-check.yaml)

**Nonparametric predictive inference** for censored and Bernoulli data:
lower and upper probabilities for future observations under Hill’s
assumption A(n), carried by M-functions.

- For data censored on the right, the left or within intervals of
  arbitrary width: the M-functions of a sharing assumption that
  reproduces the right-censored theory of Coolen and Yan (2004) exactly,
  and survival bounds for the next observation that contain the Turnbull
  estimator under every censoring pattern.
- For Bernoulli data: bounds on the number of events among several
  future observations (Coolen 1998), and the same bounds from censored
  data by completions.

Computation runs on a shared pure-C back-end that is also bound from
Python, and the two front ends are checked against each other value by
value.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mtloots/npistat", subdir = "npistat")
```

## Usage

The six-observation right-censored example of Yan (2002), in units of
1/35:

``` r
library(npistat)
tm <- c(3, 4, 6, 7, 8, 9); ev <- c(1, 0, 1, 0, 1, 1)
s <- npi_survival(tm, ifelse(ev == 1, tm, Inf), t = c(5, 7.5, 9.5))
transform(s, lower35 = lower * 35, upper35 = upper * 35)
#>     t     lower     upper lower35 upper35
#> 1 5.0 0.6857143 0.8571429      24      30
#> 2 7.5 0.4571429 0.6857143      16      24
#> 3 9.5 0.0000000 0.2285714       0       8
```

Interval-censored data, with the members of the construction and the
Turnbull set:

``` r
L <- c(3.5, 3.4, 4.7, 6.9, 5.6, 2.4); R <- c(6.5, 6.8, 5.9, 6.9, 7.9, 6.4)
npi_survival(L, R, t = 5.75, members = TRUE)
#>      t     lower upper   S.lower S.upper   P.lower P.upper   C.lower C.upper
#> 1 5.75 0.1714286     1 0.2142857       1 0.1714286       1 0.1428571       1
#>   TB.lower TB.upper
#> 1      0.2        1
```

A portfolio of 1000 obligors with two defaults: bounds on the number of
defaults among the next 100.

``` r
npi_bernoulli(n = 1000, k = 2, m = 100, r = 1:5)
#>   r        lower        upper
#> 1 1 0.1736289189 0.2488904017
#> 2 2 0.0231059533 0.0434821250
#> 3 3 0.0027297815 0.0063736954
#> 4 4 0.0003005056 0.0008384959
#> 5 5 0.0000315104 0.0001023245
```
