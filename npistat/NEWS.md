# npistat 0.1.1

* `npi_turnbull()`, and through it `npi_survival()`: the self-consistency iteration now only locates the
  support; a damped active-set Newton method on the optimality conditions then fixes the masses to machine
  precision, with one self-consistency step at the end. The iteration alone could stop short of the optimum
  by up to 1e-4 in the survival bounds (current-status and interval-only inspection data). The plain iteration
  remains the fallback when the masses are not identified. New tests check the optimality conditions to 1e-12.
* Python front: a cached build of the C back-end older than the C source is rebuilt instead of reused.

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
