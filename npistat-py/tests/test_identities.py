"""Identity tests for the Python front end.

The parity harness at the repository root proves that the R and Python fronts return byte-identical
values from the shared C. That is agreement, not correctness: both fronts could agree on a wrong
number. These tests assert the mathematics independently, mirroring
npistat/tests/testthat/test-npistat.R, so the Python package is checked on its own terms.
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import npistat as p  # noqa: E402

INF = float("inf")


def close(a, b, tol):
    return abs(a - b) <= tol * max(1.0, abs(b))


def test_six_observation_example():
    tm = [3, 4, 6, 7, 8, 9]; ev = [1, 0, 1, 0, 1, 1]
    R = [t if e else INF for t, e in zip(tm, ev)]
    s = p.npi_survival(tm, R, [1.5, 3.5, 5, 6.5, 7.5, 8.5, 9.5])
    assert [round(r["lower"] * 35, 9) for r in s] == [30, 25, 24, 18, 16, 8, 0]
    assert [round(r["upper"] * 35, 9) for r in s] == [35, 30, 30, 24, 24, 16, 8]


def test_coolen_yan_table_2():
    tA = [90, 142, 150, 269, 291, 468, 680, 837, 890, 1037, 1090, 1113, 1153, 1297, 1429, 1577]
    cA = [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1]
    tB = [272, 362, 373, 383, 519, 563, 650, 827, 919, 978, 1100, 1307, 1360, 1476]
    cB = [0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1]
    pubA = [0.0588] * 6 + [0.0053, 0.0642, 0.0642, 0.0080, 0.0722, 0.0120, 0.0168] + [0.1011] * 4
    pubB = [0.0667] * 4 + [0.0061, 0.0073, 0.0089, 0.0111, 0.1000, 0.0167, 0.0233, 0.0350, 0.1750, 0.0875, 0.2625]
    MA = p.npi_mfun(tA, [t if c == 0 else INF for t, c in zip(tA, cA)])
    MB = p.npi_mfun(tB, [t if c == 0 else INF for t, c in zip(tB, cB)])
    assert len(MA) == 17 and len(MB) == 15
    assert max(abs(m[2] - q) for m, q in zip(MA, pubA)) < 5e-5
    assert max(abs(m[2] - q) for m, q in zip(MB, pubB)) < 5e-5


def test_sharing_equals_rc_closed_form():
    rng = random.Random(1)
    for _ in range(100):
        n = rng.randint(3, 15)
        a = [rng.uniform(0, 10) for _ in range(n)]
        ev = [int(rng.random() < 0.6) for _ in range(n)]
        S = sorted(m[2] for m in p.npi_mfun(a, [x if e else INF for x, e in zip(a, ev)]))
        C = sorted(m[2] for m in p.npi_rc(a, ev))
        assert all(close(x, y, 1e-9) for x, y in zip(S, C))
        assert close(sum(S), 1.0, 1e-12)


def test_reduces_to_A_n_on_exact_data():
    x = [1.3, 2.7, 4.1, 6.6, 9.2]; n = len(x)
    for r in p.npi_survival(x, x, [0.5, 2, 5, 8, 10]):
        k = sum(v > r["t"] for v in x)
        assert close(r["lower"], k / (n + 1), 1e-12) and close(r["upper"], (k + 1) / (n + 1), 1e-12)


def test_ties_follow_H_n():
    r = p.npi_survival([1, 2, 2, 4], [1, 2, 2, 4], 3)[0]
    assert close(r["lower"], 0.2, 1e-12) and close(r["upper"], 0.4, 1e-12)


def test_containments():
    rng = random.Random(2)
    for _ in range(100):
        n = rng.randint(4, 20)
        a = [rng.uniform(0, 10) for _ in range(n)]
        ex = [rng.random() < 0.5 for _ in range(n)]
        R = [x if e else x + rng.uniform(0.5, 4) for x, e in zip(a, ex)]
        for r in p.npi_survival(a, R, [0.37 + 1.3 * j for j in range(11)], members=True):
            assert r["lower"] <= r["TB.lower"] + 1e-12 and r["TB.upper"] <= r["upper"] + 1e-12
            assert r["C.lower"] <= r["lower"] + 1e-12 and r["upper"] <= r["C.upper"] + 1e-12
            assert r["S.lower"] >= r["C.lower"] - 1e-12


def test_leak_example():
    L = [3.5, 3.4, 4.7, 6.9, 5.6, 2.4]; R = [6.5, 6.8, 5.9, 6.9, 7.9, 6.4]
    r = p.npi_survival(L, R, 5.75, members=True)[0]
    assert close(r["TB.lower"], 0.2, 1e-6) and close(r["S.lower"], 3 / 14, 1e-9)
    assert close(r["lower"], 6 / 7 * 0.2, 1e-6)
    assert [round(c[2], 6) for c in p.npi_turnbull(L, R)] == [0.8, 0.2]


def test_bernoulli():
    for n in (5, 50, 500):
        for k in range(4):
            _, lo, up = p.npi_bernoulli(n, k, 1, 1)[0]
            assert close(lo, k / (n + 1), 1e-10) and close(up, (k + 1) / (n + 1), 1e-10)
    for n in (50, 146, 394, 20000):
        assert close(p.npi_bernoulli(n, 0, 100, 1)[0][2], 100 / (n + 100), 1e-9)
    assert close(p.npi_bernoulli(2, 0, 2, 1)[0][2], 0.5, 1e-12)


def test_completion_equals_bernoulli_at_common_horizon():
    for k in range(4):
        L = [0.5] * k + [1.0] * (200 - k); R = [0.5] * k + [INF] * (200 - k)
        _, flo, fup = p.npi_future(L, R, 0.95, 50, 2, side="below")[0]
        _, blo, bup = p.npi_bernoulli(200, k, 50, 2)[0]
        assert close(flo, blo, 1e-10) and close(fup, bup, 1e-10)


def test_twolevel_closed_form_and_ceiling():
    from math import comb
    def tail(m, r, t):
        return sum(comb(m, k) * t ** k * (1 - t) ** (m - k) for k in range(r, m + 1))
    rng = random.Random(3)
    for _ in range(20):
        n = rng.randint(5, 30)
        th = [0.0 if rng.random() < 0.4 else rng.random() / 50 for _ in range(n)]
        s = sorted(th)
        for m, r in ((1, 1), (20, 1), (20, 2), (50, 5)):
            lo = sum(tail(m, r, t) for t in [0.0] + s) / (n + 1)
            up = sum(tail(m, r, t) for t in s + [1.0]) / (n + 1)
            _, a, b = p.npi_twolevel(th, th, m, r)[0]
            assert close(a, lo, 1e-10) and close(b, up, 1e-10) and b >= 1 / (n + 1) - 1e-12
