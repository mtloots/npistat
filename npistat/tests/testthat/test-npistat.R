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
    expect_true(all(s$S.lower >= s$C.lower - 1e-12))                # Proposition 4(a)
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
    for (m in c(1, 20, 500)) for (r in unique(c(1, 2, min(5, m), m))) {
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
