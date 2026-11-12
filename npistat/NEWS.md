# npistat 0.2.0

* One tie convention for every member. A censored region is read as half-open, `(L, R]`: an exact value tied
  with its right end lies inside it, an exact value tied with its left end lies outside and below it, and a
  region that ends at a value lies below a region that starts there. `npi_survival()` now computes the Turnbull
  estimator from the regions with ties broken, as it always did the sharing specification; before, the
  estimator took the regions as closed sets, so that with ties the anchored member and the sharing member
  described different data, and under right censoring with an event tied with a censoring the construction was
  not `rc-A(n)`. `npi_turnbull()` breaks ties in the same way by default (`untie = FALSE` gives the estimator
  for closed regions, as before). Bounds on data with ties can change in the third decimal.
* The right end of a censored region is now moved just past the exact values tied with it and no further.
  Before, it was moved past censored regions starting at the same value as well, so that two touching regions
  overlapped when an exact value sat at the common point and not otherwise.
* The tie-breaking increment is chosen from the data by default (`eps = 0`): the smallest gap between distinct
  endpoint values divided by `64 (n + 2)`. The bounds no longer depend on the unit of measurement; the old
  default, an absolute `1e-7`, could change the order of values closer together than `n` times that.
* `npi_mfun()` and `npi_survival()`: the masses of the sharing specification are now placed between the ends
  of the observed regions. Each mass starts at the left end of its region and ends at the right end, one unit
  starts at the origin and one ends at infinity, and left and right ends are paired in increasing order. The
  lower survival function is unchanged. The upper survival function is one unit plus the masses of the
  regions that end after `t`; it can be lower than before under interval censoring, and it never exceeds
  the completion bound, so `npi_survival()` no longer clips. Exact and right-censored data give the same
  M-functions as before. `npi_mfun()` can now return up to `2n + 1` intervals; their ends are data values, and
  an interval with `from == to` is a point mass at a tied value.
* `npi_twolevel()` and `npi_rates()` use the same construction, so their bounds can move with it.
* `npi_turnbull()`, and through it `npi_survival()`: with right-censored observations the estimator was wrong
  when a finite right endpoint exceeded three times the largest left endpoint plus one, and its upper survival
  function was zero beyond that value. The stand-in for an infinite right endpoint is now taken beyond every
  finite endpoint, and a clique that is unbounded above is returned with `right = Inf`.
* `npi_mfun()`, `npi_survival()` and `npi_turnbull()` refuse negative endpoints: the scale starts at the origin.

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
