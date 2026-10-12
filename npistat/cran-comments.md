## Submission

This is a first submission of 'npistat'.

The package computes nonparametric predictive inference (NPI): lower and upper probabilities for
future observations under Hill's assumption A(n), carried by M-functions. For censored data it gives
the M-functions of a sharing assumption that reproduces the right-censored theory of Coolen and Yan
(2004) exactly, the Turnbull (1976) estimator, and survival bounds for the next observation; for
Bernoulli data it gives bounds on the number of events among several future observations.

Points raised by the CRAN team on my earlier submissions are taken as standing requirements.

**`\value` tags.** Every `.Rd` documenting an exported function states what the function returns.
The package-level page `npistat-package.Rd` has no `\value`, which is correct for an overview page.

**References in the Description.** Hill (1968), Coolen and Yan (2004) and Coolen (1998) are cited
with their DOIs, each confirmed against Crossref rather than copied from a bibliography. Turnbull
(1976) has no DOI and is cited by name.

Software and package names are single-quoted in the Description ('Python', 'npistat').

## Test environments

	local macOS (arm64), R 4.5, `R CMD check --as-cran` on the built tarball

Win-builder results will be recorded here once they have been run on this tarball; nothing is
claimed before it has run.

## R CMD check results

0 errors | 0 warnings | 2 notes locally

The first note is "New submission". Until the GitHub repository is public it also lists the
repository and issues URLs as not found; they will resolve before submission.

The second note is about the local 'tidy' and 'V8' installations, a property of the maintainer's
machine rather than of the package.

## Notes for the reviewer

* The package contains compiled C. It uses only libm and R's own headers, needs no system
  requirement beyond a C compiler, and registers its native routines with `R_registerRoutines`
  and `R_useDynamicSymbols(dll, FALSE)`.
* The same C sources back a 'Python' package of the same name. The two front ends are compared
  value by value by a parity harness in the repository, and the tests here encode published values
  and proved properties (the 32 M-function values of Coolen and Yan's Table 2, Hill's bounds,
  containment of the Turnbull estimator) rather than recording current output.
* No example, test or vignette uses more than one core.
