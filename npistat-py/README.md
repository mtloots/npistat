# npistat

Nonparametric predictive inference (NPI) for censored and Bernoulli data, from one pure-C back-end
shared by an R package and a Python package of the same name.

* M-functions and survival bounds for the next observation from data censored on the right, the left
  or within intervals of arbitrary width. Under right censoring the M-functions are the rc-A(n) of
  Coolen and Yan (2004), reproduced exactly; under every censoring pattern the bounds contain the
  Turnbull estimator.
* Bounds on the number of events among several future Bernoulli observations (Coolen 1998), and the
  same bounds from censored data by completions.

```r
tm <- c(3, 4, 6, 7, 8, 9); ev <- c(1, 0, 1, 0, 1, 1)
npistat::npi_survival(tm, ifelse(ev == 1, tm, Inf), t = c(5, 7.5, 9.5))
npistat::npi_bernoulli(n = 1000, k = 2, m = 100, r = 1:5)
```

```python
import npistat
npistat.npi_survival([3, 4, 6, 7, 8, 9], [3, float("inf"), 6, float("inf"), 8, 9], [5, 7.5, 9.5])
npistat.npi_bernoulli(1000, 2, 100, [1, 2, 3])
```
