# npistat

[![R-CMD-check](https://github.com/mtloots/npistat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/npistat/actions/workflows/R-CMD-check.yaml)

Nonparametric predictive inference (NPI) for censored and Bernoulli data: lower and upper
probabilities for future observations under Hill's assumption A(n), carried by M-functions. For data
censored on the right, the left or within intervals of arbitrary width, a sharing assumption that
reproduces the right-censored theory of Coolen and Yan (2004) exactly, and survival bounds that
contain the Turnbull estimator under every censoring pattern; for Bernoulli data, bounds on the number
of events among several future observations.

One shared pure-C back-end (`npistat.c`) serves two front-ends:

	npistat/      R package    (.C bindings)
	npistat-py/   Python package (ctypes; pip builds the C at install time)

	parity_npistat.sh      every exported routine, byte-identical between the two fronts
	parity_sync_guard.sh   refuses a parity run against a stale copy of the C

Companion software to the manuscripts "Nonparametric predictive inference for interval-censored data:
M-functions anchored to the self-consistent estimator" and "Predictive bounds for low-default
portfolios" (M. T. Loots, in preparation).
