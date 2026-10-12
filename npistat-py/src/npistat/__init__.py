"""npistat: nonparametric predictive inference for censored and Bernoulli data.

Lower and upper probabilities for future observations under Hill's assumption A(n): M-functions and
survival bounds for data censored on the right, the left or within intervals (npi_mfun,
npi_survival, npi_rc, npi_turnbull), and bounds for events among several future Bernoulli
observations (npi_bernoulli, npi_future). ctypes binding to a shared C back-end also used by the R
package 'npistat'. Pure standard library.
"""
from ._core import (npi_bernoulli, npi_future, npi_mfun, npi_rc, npi_turnbull, npi_survival)
__all__ = ["npi_bernoulli", "npi_future", "npi_mfun", "npi_rc", "npi_turnbull", "npi_survival"]
__version__ = "0.1.0"
