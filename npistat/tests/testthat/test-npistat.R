## The tests check published values and proved properties, not merely that the code runs.

test_that("the six-observation example of Yan (2002) is reproduced in units of 1/35", {
  tm <- c(3, 4, 6, 7, 8, 9); ev <- c(1, 0, 1, 0, 1, 1)
  s <- npi_survival(tm, ifelse(ev == 1, tm, Inf), c(1.5, 3.5, 5, 6.5, 7.5, 8.5, 9.5))
  expect_equal(s$lower * 35, c(30, 25, 24, 18, 16, 8, 0), tolerance = 1e-9)
  expect_equal(s$upper * 35, c(35, 30, 30, 24, 24, 16, 8), tolerance = 1e-9)
})

test_that("the 32 published M-function values of Coolen and Yan (2004, Table 2) are reproduced", {
  tA <- c(90,142,150,269,291,468,680,837,890,1037,1090,1113,1153,1297,1429,1577)
  cA <- c(0,0,0,0,0,1,0,0,1,0,1,1,0,0,0,1)
  tB <- c(272,362,373,383,519,563,650,827,919,978,1100,1307,1360,1476)
  cB <- c(0,0,0,1,1,1,1,0,1,1,1,0,1,1)
  pubA <- c(rep(0.0588,6),0.0053,0.0642,0.0642,0.0080,0.0722,0.0120,0.0168,rep(0.1011,4))
  pubB <- c(rep(0.0667,4),0.0061,0.0073,0.0089,0.0111,0.1000,0.0167,0.0233,0.0350,0.1750,0.0875,0.2625)
  MA <- npi_mfun(tA, ifelse(cA == 0, tA, Inf)); MB <- npi_mfun(tB, ifelse(cB == 0, tB, Inf))
  expect_equal(nrow(MA), 17); expect_equal(nrow(MB), 15)
  expect_lt(max(abs(MA$mass - pubA)), 5e-5); expect_lt(max(abs(MB$mass - pubB)), 5e-5)
})

test_that("the sharing M-functions equal the rc-A(n) closed form on right-censored data", {
  set.seed(1)
  for (rep in 1:100) {
    n <- sample(3:15, 1); a <- runif(n, 0, 10); ev <- as.integer(runif(n) < 0.6)
    S <- npi_mfun(a, ifelse(ev == 1, a, Inf)); C <- npi_rc(a, ev)
    expect_equal(sort(S$mass), sort(C$mass), tolerance = 1e-9)
    expect_equal(sum(S$mass), 1, tolerance = 1e-12)
  }
})

test_that("the bounds reduce to A(n) on exact data", {
  x <- c(1.3, 2.7, 4.1, 6.6, 9.2); n <- length(x); t <- c(0.5, 2, 5, 8, 10)
  s <- npi_survival(x, x, t)
  expect_equal(s$lower, sapply(t, function(u) sum(x > u)) / (n + 1), tolerance = 1e-12)
  expect_equal(s$upper, sapply(t, function(u) sum(x > u) + 1) / (n + 1), tolerance = 1e-12)
})

test_that("ties follow Hill's H(n)", {
  s <- npi_survival(c(1, 2, 2, 4), c(1, 2, 2, 4), 3)
  expect_equal(c(s$lower, s$upper), c(1, 2) / 5, tolerance = 1e-12)
})

test_that("the construction contains the Turnbull set and lies inside the completion bound", {
  set.seed(2)
  for (rep in 1:100) {
    n <- sample(4:20, 1); a <- runif(n, 0, 10); ex <- runif(n) < 0.5
    R <- ifelse(ex, a, a + runif(n, 0.5, 4)); t <- seq(0.37, 13.37, by = 1.3)
    s <- npi_survival(a, R, t, members = TRUE)
    expect_true(all(s$lower <= s$TB.lower + 1e-12 & s$TB.upper <= s$upper + 1e-12))
    expect_true(all(s$C.lower <= s$lower + 1e-12 & s$upper <= s$C.upper + 1e-12))
    expect_true(all(s$S.lower >= s$C.lower - 1e-12 & s$S.upper <= s$C.upper + 1e-12))  # the sharing member inside M_C
    expect_true(all(s$P.lower >= s$C.lower - 1e-12 & s$P.upper <= s$C.upper + 1e-12))  # and the anchored member
  }
})

test_that("the Turnbull estimator satisfies its optimality conditions to machine precision", {
  gen <- function(n, pex) { Tt <- rexp(n, 1/5); ex <- runif(n) < pex; L <- Tt; R <- Tt
    for (i in which(!ex)) { ins <- sort(runif(sample(2:6, 1), 0, 15)); L[i] <- max(c(0, ins[ins < Tt[i]])); hi <- ins[ins >= Tt[i]]; R[i] <- if (length(hi)) min(hi) else Inf }
    list(L = L, R = R) }
  set.seed(7)
  for (rep in 1:30) {                       # inspection designs, including interval-only data at n = 240, where the
    n <- c(10, 40, 240)[1 + rep %% 3]; d <- gen(n, if (rep %% 2) 0.5 else 0)    # plain iteration stopped up to 1e-4 short
    tb <- npi_turnbull(d$L, d$R)                       # an unbounded clique has right = Inf, inside a right-censored region
    A <- outer(seq_len(n), seq_len(nrow(tb)), function(i, j) tb$left[j] >= d$L[i] & tb$right[j] <= d$R[i]) * 1
    g <- as.vector(t(A) %*% (1 / as.vector(A %*% tb$mass))); on <- tb$mass > 0
    expect_lt(max(abs(g[on] - n)) / n, 1e-12)                      # sum_i a_ij / d_i = n on the support
    if (any(!on)) expect_lte(max(g[!on] - n) / n, 1e-12)            # and at most n off it
    expect_equal(sum(tb$mass), 1, tolerance = 1e-14)
  }
})

test_that("the leak example: sharing exceeds the Turnbull lower bound, the construction does not", {
  L <- c(3.5, 3.4, 4.7, 6.9, 5.6, 2.4); R <- c(6.5, 6.8, 5.9, 6.9, 7.9, 6.4)
  s <- npi_survival(L, R, 5.75, members = TRUE)
  expect_equal(s$TB.lower, 0.2, tolerance = 1e-6)
  expect_equal(s$S.lower, 3/14, tolerance = 1e-9)
  expect_equal(s$lower, 6/7 * 0.2, tolerance = 1e-6)
  tb <- npi_turnbull(L, R)
  expect_equal(tb$mass, c(0.8, 0.2), tolerance = 1e-6)
})

test_that("Bernoulli bounds: Hill at m = 1, m/(n+m) at zero events, counting identity", {
  for (n in c(5, 50, 500)) for (k in 0:3) {
    b <- npi_bernoulli(n, k, 1, 1)
    expect_equal(c(b$lower, b$upper), c(k, k + 1) / (n + 1), tolerance = 1e-10)
  }
  for (n in c(50, 146, 394, 20000)) expect_equal(npi_bernoulli(n, 0, 100, 1)$upper, 100 / (n + 100), tolerance = 1e-9)
  expect_equal(npi_bernoulli(2, 0, 2, 1)$upper, 1 - 2/4, tolerance = 1e-12)   # 2/3 x 3/4 = 2/4 survive
  b <- npi_bernoulli(1000, 2, 100, 1:10)
  expect_true(all(diff(b$lower) <= 1e-15) && all(diff(b$upper) <= 1e-15) && all(b$lower <= b$upper))
})

test_that("completion bounds for a portfolio censored at a common horizon equal the Bernoulli bounds", {
  for (k in 0:3) {
    L <- c(rep(0.5, k), rep(1, 200 - k)); R <- c(rep(0.5, k), rep(Inf, 200 - k))
    f <- npi_future(L, R, 0.95, m = 50, r = 2, side = "below"); b <- npi_bernoulli(200, k, 50, 2)
    expect_equal(c(f$lower, f$upper), c(b$lower, b$upper), tolerance = 1e-10)
  }
})

test_that("two-level bounds on known rates equal the closed form, and the 1/(n+1) ceiling holds", {
  closed <- function(th, m, r) { s <- sort(th)
    c(mean(pbinom(r - 1, m, c(0, s), lower.tail = FALSE)), mean(pbinom(r - 1, m, c(s, 1), lower.tail = FALSE))) }
  set.seed(3)
  for (rep in 1:30) { n <- sample(5:40, 1); th <- ifelse(runif(n) < 0.4, 0, runif(n) / 50)
    for (m in c(1, 20, 500)) for (r in unique(pmin(c(1, 2, 5, m), m))) {
      b <- npi_twolevel(th, th, m, r)
      expect_equal(c(b$lower, b$upper), closed(th, m, r), tolerance = 1e-10)
      expect_gte(b$upper, 1 / (n + 1) - 1e-12) } }
})

test_that("two-level with one unit is A(n) on the rates' expectation; interval years widen nothing below zero", {
  th <- c(0, 0.01, 0.02, 0.05); b <- npi_twolevel(th, th, 1, 1)
  expect_equal(c(b$lower, b$upper), c(sum(th) / 5, (sum(th) + 1) / 5), tolerance = 1e-12)
  iv <- npi_rates(c(0, 1, 2), c(1000, 1000, 1000))
  expect_equal(iv$L, c(0, 1, 2) / 1001); expect_equal(iv$R, c(1, 2, 3) / 1001)
  b <- npi_twolevel(iv$L, iv$R, 100, 1:3)
  expect_true(all(b$lower >= 0 & b$lower <= b$upper & b$upper <= 1))
})

test_that("the sharing M-function has proper intervals, sums to one, and is A(n) on exact data", {
  x <- c(1.3, 2.7, 4.1, 6.6, 9.2); M <- npi_mfun(x, x)
  expect_equal(M$from, c(0, x), tolerance = 1e-6); expect_equal(M$to, c(x, Inf), tolerance = 1e-6); expect_equal(M$mass, rep(1 / 6, 6), tolerance = 1e-12)
  set.seed(5)
  for (rep in 1:200) {
    n <- sample(2:15, 1); a <- round(runif(n, 1, 10), 2); ty <- sample(c("e", "i", "r", "l"), n, TRUE, prob = c(.35, .4, .15, .1))
    L <- ifelse(ty == "l", 0, a); R <- ifelse(ty == "e", a, ifelse(ty == "r", Inf, round(a + runif(n, 0.3, 7), 2)))
    if (anyDuplicated(c(a, R[ty == "i"]))) next
    M <- npi_mfun(L, R)
    expect_equal(sum(M$mass), 1, tolerance = 1e-12); expect_true(all(M$from < M$to)); expect_lte(nrow(M), 2 * n + 1)
  }
})

test_that("on nested regions the sharing M-function is the ordering rule: every consistent ordering equally likely", {
  perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(k) lapply(perms(v[-k]), function(p) c(v[k], p))))
  ## enumeration: a consistent ordering has every left endpoint below every later right endpoint; each of its n + 1
  ## gaps lies between the largest left endpoint below it and the smallest right endpoint above it
  rule <- function(L, R) { n <- length(L); out <- NULL; E <- 0
    for (p in perms(seq_len(n))) { ok <- TRUE
      for (k in seq_len(n - 1)) if (max(L[p[1:k]]) >= min(R[p[(k + 1):n]])) { ok <- FALSE; break }
      if (!ok) next
      E <- E + 1
      for (k in 0:n) out <- rbind(out, c(if (k == 0) 0 else max(L[p[1:k]]), if (k == n) Inf else min(R[p[(k + 1):n]]))) }
    d <- data.frame(from = out[, 1], to = out[, 2], mass = 1 / ((n + 1) * E)); d <- stats::aggregate(mass ~ from + to, d, sum)
    d[order(d$from, d$to), ] }
  nested <- function(L, R) { n <- length(L)
    for (a in 1:n) for (b in 1:n) if (a < b && L[a] < R[a] && L[b] < R[b] &&
      !((L[a] <= L[b] && R[b] <= R[a]) || (L[b] <= L[a] && R[a] <= R[b]) || R[a] < L[b] || R[b] < L[a])) return(FALSE)
    TRUE }
  set.seed(8); done <- 0
  while (done < 60) { n <- sample(2:6, 1); a <- round(runif(n, 1, 10), 2); ty <- sample(c("e", "i", "r"), n, TRUE, prob = c(.4, .45, .15))
    L <- a; R <- ifelse(ty == "e", a, ifelse(ty == "r", Inf, round(a + runif(n, 0.3, 7), 2)))
    if (anyDuplicated(c(L, R[ty == "i"])) || !nested(L, R)) next
    done <- done + 1; A <- npi_mfun(L, R); B <- rule(L, R)
    expect_equal(nrow(A), nrow(B)); expect_equal(A$from, B$from, tolerance = 1e-5); expect_equal(A$to, B$to, tolerance = 1e-5)
    expect_equal(A$mass, B$mass, tolerance = 1e-12)
    t <- sort(unique(c(L, R[is.finite(R)]))); t <- c(t[1] / 2, (t[-1] + t[-length(t)]) / 2, max(t) + 1); s <- npi_survival(L, R, t, members = TRUE)
    expect_equal(s$S.lower, sapply(t, function(u) sum(B$mass[B$from >= u])), tolerance = 1e-12)
    expect_equal(s$S.upper, sapply(t, function(u) sum(B$mass[B$to > u])), tolerance = 1e-12)
  }
})

test_that("right-censored observations: the Turnbull estimator when a finite right endpoint lies far beyond the left endpoints, and beyond all the data", {
  ## an exact value 0.2, an interval [0.5, 20] and a right-censored observation at 1: cliques {0.2} and [1, 20]
  tb <- npi_turnbull(c(0.2, 0.5, 1), c(0.2, 20, Inf))
  expect_equal(tb$left, c(0.2, 1)); expect_equal(tb$right, c(0.2, 20)); expect_equal(tb$mass, c(1, 2) / 3, tolerance = 1e-12)
  s <- npi_survival(c(0.2, 0.5, 1), c(0.2, 20, Inf), c(10, 30), members = TRUE)
  expect_equal(s$TB.upper, c(2 / 3, 0), tolerance = 1e-12); expect_equal(s$TB.lower, c(0, 0), tolerance = 1e-12)
  ## an event at 2 and a censoring at 3: the Kaplan-Meier estimator keeps one half beyond 3, however far out it is read
  tb <- npi_turnbull(c(2, 3), c(2, Inf)); expect_equal(tb$right, c(2, Inf)); expect_equal(tb$mass, c(0.5, 0.5), tolerance = 1e-12)
  s <- npi_survival(c(2, 3), c(2, Inf), c(2.5, 100, 1e6), members = TRUE)
  expect_equal(s$TB.upper, rep(0.5, 3), tolerance = 1e-12); expect_equal(s$P.upper, rep(2 / 3, 3), tolerance = 1e-12)
  expect_equal(s$upper, rep(2 / 3, 3), tolerance = 1e-12)
  ## random right-censored data: the Turnbull upper survival function is the Kaplan-Meier estimator at every point, far out included
  set.seed(12)
  for (rep in 1:50) { n <- sample(3:15, 1); tm <- runif(n, 0, 10); ev <- runif(n) < 0.6; t <- c(sort(tm)[-1] - diff(sort(tm)) / 2, max(tm) + c(1, 50, 1e4))
    km <- sapply(t, function(u) { o <- order(tm); s <- 1; for (i in o) if (ev[i] && tm[i] <= u) s <- s * (1 - 1 / sum(tm >= tm[i])); s })
    expect_equal(npi_survival(tm, ifelse(ev, tm, Inf), t, members = TRUE)$TB.upper, km, tolerance = 1e-10) }
})

test_that("ties: one convention for every member, the half-open reading of a censored region", {
  ## events at 1, 1, 1, right-censorings at 1, 1, an event at 2: the estimator is Kaplan-Meier (an event precedes a
  ## censoring tied with it) and the construction is rc-A(n), whose lower survival probability at 1.5 is 2/7
  L <- c(1, 1, 1, 1, 1, 2); R <- c(1, 1, 1, Inf, Inf, 2)
  tb <- npi_turnbull(L, R); expect_equal(tb$left, c(1, 2)); expect_equal(tb$mass, c(0.5, 0.5), tolerance = 1e-12)
  s <- npi_survival(L, R, 1.5, members = TRUE)
  expect_equal(s$lower, 2 / 7, tolerance = 1e-12); expect_equal(s$S.lower, 2 / 7, tolerance = 1e-12); expect_equal(s$TB.upper, 0.5, tolerance = 1e-12)
  expect_equal(npi_turnbull(L, R, untie = FALSE)$mass, c(0.75, 0.25), tolerance = 1e-10)       # closed regions, as given
  ## tied right-censored data in general: Turnbull is Kaplan-Meier, and the anchored member lies inside rc-A(n)
  set.seed(21)
  for (rep in 1:300) { n <- sample(2:8, 1); tm <- sample(1:4, n, TRUE); ev <- runif(n) < 0.6; t <- c(0.5, 1.5, 2.5, 3.5, 4.5, 50)
    km <- sapply(t, function(u) { s <- 1; for (v in sort(unique(tm[ev]))) if (v <= u) s <- s * (1 - sum(ev & tm == v) / sum(tm >= v)); s })
    s <- npi_survival(tm, ifelse(ev, tm, Inf), t, members = TRUE)
    expect_equal(s$TB.upper, km, tolerance = 1e-10)
    expect_true(all(s$P.lower >= s$S.lower - 1e-12 & s$P.upper <= s$S.upper + 1e-12))
    expect_equal(s$lower, s$S.lower, tolerance = 1e-12); expect_equal(s$upper, s$S.upper, tolerance = 1e-12) }
  ## an exact value at the end of an interval lies inside it, at the start outside and below it
  expect_equal(npi_turnbull(c(2, 5), c(5, 5))$mass, 1); expect_equal(npi_turnbull(c(2, 5), c(5, 5))$left, 5)
  tb <- npi_turnbull(c(2, 2), c(2, 5)); expect_equal(tb$left, c(2, 2)); expect_equal(tb$right, c(2, 5)); expect_equal(tb$mass, c(0.5, 0.5), tolerance = 1e-12)
  tb <- npi_turnbull(c(2, 5), c(5, 5), closedR = FALSE); expect_equal(tb$right, c(5, 5)); expect_equal(tb$mass, c(0.5, 0.5), tolerance = 1e-12)
  ## a region that ends at a value lies below a region that starts there, with or without an exact value at that point
  tb <- npi_turnbull(c(1, 3), c(3, 5)); expect_equal(tb$left, c(1, 3)); expect_equal(tb$right, c(3, 5)); expect_equal(tb$mass, c(0.5, 0.5), tolerance = 1e-12)
  tb <- npi_turnbull(c(1, 3, 3), c(3, 5, 3)); expect_equal(tb$left, c(3, 3)); expect_equal(tb$right, c(3, 5)); expect_equal(tb$mass, c(2, 1) / 3, tolerance = 1e-12)
})

test_that("the bounds do not depend on the unit of measurement, the ends of the M-function are data values, negative values are refused", {
  set.seed(22)
  for (rep in 1:100) { n <- sample(2:12, 1); a <- sample(1:8, n, TRUE); ty <- sample(c("e", "i", "r", "l"), n, TRUE, prob = c(.35, .4, .15, .1))
    L <- ifelse(ty == "l", 0, a); R <- ifelse(ty == "e", a, ifelse(ty == "r", Inf, a + sample(1:4, n, TRUE))); t <- seq(0.5, 13.5, by = 1)
    s1 <- npi_survival(L, R, t, members = TRUE)
    for (k in c(1e-6, 1e6)) expect_equal(as.matrix(npi_survival(L * k, R * k, t * k, members = TRUE)[, -1]), as.matrix(s1[, -1]), tolerance = 1e-12)
    M <- npi_mfun(L, R); expect_true(all(M$from %in% c(0, L, R)) && all(M$to %in% c(L, R, Inf))); expect_equal(sum(M$mass), 1, tolerance = 1e-12)
    expect_true(all(s1$C.lower <= s1$lower + 1e-12 & s1$upper <= s1$C.upper + 1e-12 & s1$lower <= s1$TB.lower + 1e-12 & s1$TB.upper <= s1$upper + 1e-12)) }
  expect_error(npi_survival(c(-2, -1, 3), c(-2, -1, 3), 0), "origin"); expect_error(npi_mfun(c(-1, 2), c(1, 2)), "origin")
})
