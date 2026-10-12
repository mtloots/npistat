# npistat 0.1.0

First release.

* `npi_mfun()`, `npi_survival()`: M-functions and survival bounds for data censored on the right,
  the left or within intervals; the sharing specification reproduces the rc-A(n) M-functions of
  Coolen and Yan (2004) exactly, and the survival bounds contain the Turnbull estimator.
* `npi_rc()`: the rc-A(n) M-functions in closed form.
* `npi_turnbull()`: the Turnbull estimator by self-consistency.
* `npi_bernoulli()`, `npi_future()`: bounds on the number of events among several future
  observations, from Bernoulli data and from censored data by completions.
* `npi_twolevel()`, `npi_rates()`: two-level bounds, with Hill's A(n) over exchangeable periods
  (carrying variation between periods, such as systematic default risk) and binomial counts within a
  period; a period's rate may be an interval derived from its counts.
