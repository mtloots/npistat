## npistat: R front to the shared C back-end. Every computation is in src/npistat.c; this file
## only validates input, calls .C and shapes the output.

npi_regions <- function(L, R) {
  L <- as.double(L); R <- as.double(R)
  if (length(L) != length(R)) stop("L and R must have the same length")
  if (any(is.na(L)) || any(is.na(R))) stop("L and R must not contain NA")
  if (any(R < L)) stop("every region needs L <= R")
  if (any(!is.finite(L))) stop("L must be finite; code a left-censored observation by L = 0")
  list(L = L, R = R, ex = as.integer(L == R), n = length(L))
}

#' Bounds for events among future Bernoulli observations
#'
#' Nonparametric predictive lower and upper probabilities that at least \code{r} of \code{m}
#' future observations are events, given \code{k} events among \code{n} exchangeable past
#' observations, under Hill's assumptions \eqn{A_{(n)}, \dots, A_{(n+m-1)}} and the latent-variable
#' representation of Coolen (1998). For \code{m = 1} the bounds are \eqn{k/(n+1)} and
#' \eqn{(k+1)/(n+1)}; with \code{k = 0} the upper probability of at least one event is
#' \eqn{m/(n+m)}.
#'
#' The future observations are exchangeable with the past ones and with each other, and are
#' dependent given the data; no common factor, such as systematic default risk, is represented.
#'
#' @param n number of past observations.
#' @param k number of events among them.
#' @param m number of future observations.
#' @param r vector of event counts; the bounds are for at least \code{r} events.
#' @return A data frame with columns \code{r}, \code{lower} and \code{upper}.
#' @references Coolen, F. P. A. (1998). Low structure imprecise predictive inference for Bayes'
#'   problem. Statistics and Probability Letters 36, 349--357.
#' @examples
#' npi_bernoulli(n = 1000, k = 2, m = 100, r = 1:5)
#' npi_bernoulli(n = 50, k = 0, m = 100, r = 1)   # upper = 100/150
#' @export
npi_bernoulli <- function(n, k, m, r = 1) {
  n <- as.integer(n); k <- as.integer(k); m <- as.integer(m); r <- as.integer(r)
  if (length(n) != 1 || length(k) != 1 || length(m) != 1) stop("n, k and m must be single integers")
  if (n < 1 || k < 0 || k > n || m < 1) stop("need n >= 1, 0 <= k <= n and m >= 1")
  z <- .C("C_npi_bern", n, k, m, r, length(r), lo = double(length(r)), up = double(length(r)),
          NAOK = TRUE, PACKAGE = "npistat")
  data.frame(r = r, lower = z$lo, upper = z$up)
}

#' Bounds for m future observations from censored data, by completions
#'
#' Lower and upper probabilities that at least \code{r} of \code{m} future observations fall
#' above \code{t} (\code{side = "above"}) or at or below it (\code{side = "below"}), from
#' observations known only to lie in regions \eqn{[L_i, R_i]}, taking the extremes of the
#' allocation count of Coolen (1998) over all completions of the regions to exact values.
#' Sharp when no region straddles \code{t}, as for a portfolio censored at a common horizon;
#' wide otherwise.
#'
#' @param L,R region endpoints; \code{R = Inf} for right censoring, \code{L == R} for an exact value.
#' @param t vector of thresholds.
#' @param m number of future observations.
#' @param r count; the bounds are for at least \code{r} of the \code{m}.
#' @param side \code{"above"} or \code{"below"}.
#' @return A data frame with columns \code{t}, \code{lower} and \code{upper}.
#' @examples
#' ## 1000 obligors observed to a horizon of 1, two defaults inside it
#' L <- c(0.5, 0.5, rep(1, 998)); R <- c(0.5, 0.5, rep(Inf, 998))
#' npi_future(L, R, t = 0.99, m = 100, r = 1, side = "below")
#' @export
npi_future <- function(L, R, t, m, r = 1, side = c("above", "below")) {
  side <- match.arg(side); g <- npi_regions(L, R); t <- as.double(t)
  z <- .C("C_npi_mfuture", g$n, g$L, g$R, as.integer(m), as.integer(r),
          as.integer(side == "below"), length(t), t, lo = double(length(t)), up = double(length(t)),
          NAOK = TRUE, PACKAGE = "npistat")
  data.frame(t = t, lower = z$lo, upper = z$up)
}

#' M-functions of the sharing assumption
#'
#' The M-functions of the sharing specification \eqn{M_S} for observations known only to lie in
#' regions \eqn{[L_i, R_i]}: a censored observation shares its predictive mass equally among the
#' observations whose region lies inside its own, itself included, in decreasing order of the size
#' of that set. Under right censoring these are the \eqn{\text{rc-}A_{(n)}} M-functions of Coolen
#' and Yan (2004). Ties are broken by Hill's \eqn{H_{(n)}} convention.
#'
#' @param L,R region endpoints; \code{R = Inf} for right censoring, \code{L == R} for an exact value.
#' @param closedR logical; extend the right endpoint of a censored region past an exact value tied
#'   with it (the default), or leave it open.
#' @param eps tie-breaking increment.
#' @return A data frame of intervals \code{from}, \code{to} and their \code{mass}, merged where an
#'   interval occurs twice; the masses sum to one.
#' @references Coolen, F. P. A. and Yan, K. J. (2004). Nonparametric predictive inference with
#'   right-censored data. Journal of Statistical Planning and Inference 126, 25--54.
#' @examples
#' tm <- c(3, 4, 6, 7, 8, 9); ev <- c(1, 0, 1, 0, 1, 1)
#' M <- npi_mfun(tm, ifelse(ev == 1, tm, Inf))
#' M$mass * 35
#' @export
npi_mfun <- function(L, R, closedR = TRUE, eps = 1e-7) {
  g <- npi_regions(L, R)
  u <- .C("C_npi_untie", g$n, g$L, g$R, g$ex, as.double(eps), as.integer(closedR), 1e-9,
          L2 = double(g$n), R2 = double(g$n), NAOK = TRUE, PACKAGE = "npistat")
  own <- .C("C_npi_share", g$n, u$L2, u$R2, g$ex, own = double(g$n + 1), NAOK = TRUE, PACKAGE = "npistat")$own
  p <- .C("C_npi_place", g$n, u$L2, u$R2, g$ex, own, lo = double(g$n + 1), hi = double(g$n + 1),
          m = double(g$n + 1), NAOK = TRUE, PACKAGE = "npistat")
  d <- data.frame(from = p$lo, to = p$hi, mass = p$m)
  d$from <- round(d$from, 6); d$to <- ifelse(is.finite(d$to), round(d$to, 6), Inf)
  d <- stats::aggregate(mass ~ from + to, d, sum)
  d[order(d$from, d$to), , drop = FALSE]
}

#' The right-censored M-functions in closed form
#'
#' The \eqn{\text{rc-}A_{(n)}} M-functions of Coolen and Yan (2004) for right-censored data without
#' ties, from their closed form.
#'
#' @param time observed times.
#' @param status 1 for an event, 0 for a right censoring.
#' @return A data frame of intervals \code{from}, \code{to} and their \code{mass}, in time order.
#' @examples
#' npi_rc(c(3, 4, 6, 7, 8, 9), c(1, 0, 1, 0, 1, 1))
#' @export
npi_rc <- function(time, status) {
  time <- as.double(time); status <- as.integer(status); n <- length(time)
  if (length(status) != n) stop("time and status must have the same length")
  if (anyDuplicated(time)) stop("the closed form assumes no ties; use npi_mfun()")
  z <- .C("C_npi_rc", n, time, status, lo = double(n + 1), hi = double(n + 1), m = double(n + 1),
          NAOK = TRUE, PACKAGE = "npistat")
  data.frame(from = z$lo, to = z$hi, mass = z$m)
}

#' The Turnbull estimator
#'
#' The nonparametric maximum likelihood estimator for censored data (Turnbull 1976): its maximal
#' cliques and their masses, by the self-consistency iteration.
#'
#' @param L,R region endpoints; \code{R = Inf} for right censoring.
#' @param maxit,tol iteration limit and convergence tolerance.
#' @return A data frame of cliques \code{left}, \code{right} and \code{mass}.
#' @references Turnbull, B. W. (1976). The empirical distribution function with arbitrarily
#'   grouped, censored and truncated data. JRSS B 38, 290--295.
#' @examples
#' npi_turnbull(c(3.5, 3.4, 4.7, 6.9, 5.6, 2.4), c(6.5, 6.8, 5.9, 6.9, 7.9, 6.4))
#' @export
npi_turnbull <- function(L, R, maxit = 20000L, tol = 1e-13) {
  g <- npi_regions(L, R)
  z <- .C("C_npi_turnbull", g$n, g$L, g$R, as.integer(maxit), as.double(tol), K = integer(1),
          ql = double(g$n), qr = double(g$n), p = double(g$n), NAOK = TRUE, PACKAGE = "npistat")
  k <- seq_len(z$K)
  data.frame(left = z$ql[k], right = z$qr[k], mass = z$p[k])
}

#' Survival bounds for the next observation from censored data
#'
#' Lower and upper survival probabilities \eqn{P(X_{n+1} > t)} from observations known only to lie
#' in regions \eqn{[L_i, R_i]}, by the construction of Loots: the union of the sharing
#' specification \eqn{M_S} and the specification \eqn{M_P} anchored to the Turnbull estimator,
#' clipped to the completion bound \eqn{M_C}. On exact data the bounds are those of \eqn{A_{(n)}},
#' on right-censored data those of \eqn{\text{rc-}A_{(n)}}, and they contain the Turnbull estimator
#' under every censoring pattern. Evaluate away from the data values.
#'
#' @param L,R region endpoints; \code{R = Inf} for right censoring, \code{L == R} for an exact value.
#' @param t vector of evaluation points.
#' @param closedR,eps tie convention, as in \code{\link{npi_mfun}}.
#' @param members logical; also return the bounds of the three members and the Turnbull set.
#' @return A data frame with \code{t}, \code{lower}, \code{upper} and, if \code{members}, the
#'   columns \code{S.lower}, \code{S.upper}, \code{P.lower}, \code{P.upper}, \code{C.lower},
#'   \code{C.upper}, \code{TB.lower} and \code{TB.upper}.
#' @examples
#' L <- c(3.5, 3.4, 4.7, 6.9, 5.6, 2.4); R <- c(6.5, 6.8, 5.9, 6.9, 7.9, 6.4)
#' npi_survival(L, R, t = c(5.75, 7.5), members = TRUE)
#' @export
npi_survival <- function(L, R, t, closedR = TRUE, eps = 1e-7, members = FALSE) {
  g <- npi_regions(L, R); t <- as.double(t); nt <- length(t)
  z <- .C("C_npi_bounds", g$n, g$L, g$R, g$ex, as.double(eps), as.integer(closedR), nt, t,
          lo = double(nt), up = double(nt), slo = double(nt), sup = double(nt), plo = double(nt),
          pup = double(nt), clo = double(nt), cup = double(nt), tlo = double(nt), tup = double(nt),
          NAOK = TRUE, PACKAGE = "npistat")
  out <- data.frame(t = t, lower = z$lo, upper = z$up)
  if (members) out <- cbind(out, S.lower = z$slo, S.upper = z$sup, P.lower = z$plo, P.upper = z$pup,
                            C.lower = z$clo, C.upper = z$cup, TB.lower = z$tlo, TB.upper = z$tup)
  out
}

#' npistat: nonparametric predictive inference for censored and Bernoulli data
#'
#' Lower and upper probabilities for future observations under Hill's assumption
#' \eqn{A_{(n)}}, from one pure-C back-end shared with the Python package of the same name:
#' M-functions and survival bounds for data censored on the right, the left or within intervals
#' (\code{\link{npi_mfun}}, \code{\link{npi_survival}}, \code{\link{npi_rc}},
#' \code{\link{npi_turnbull}}), and bounds for events among several future Bernoulli observations
#' (\code{\link{npi_bernoulli}}, \code{\link{npi_future}}).
#' @useDynLib npistat, .registration = TRUE
#' @keywords internal
"_PACKAGE"
