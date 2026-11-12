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


def test_turnbull_optimality_to_machine_precision():
    import random
    rng = random.Random(7)
    def gen(n, pex):
        L, R = [], []
        for _ in range(n):
            t = rng.expovariate(1 / 5)
            if rng.random() < pex: L.append(t); R.append(t); continue
            ins = sorted(rng.uniform(0, 15) for _ in range(rng.randint(2, 6)))
            lo = [x for x in ins if x < t]; hi = [x for x in ins if x >= t]
            L.append(max(lo) if lo else 0.0); R.append(min(hi) if hi else INF)
        return L, R
    for rep in range(30):
        n = (10, 40, 240)[rep % 3]; L, R = gen(n, 0.5 if rep % 2 else 0.0)
        tb = p.npi_turnbull(L, R)                     # an unbounded clique has right = inf, inside a right-censored region
        A = [[1.0 if (c[0] >= L[i] and c[1] <= R[i]) else 0.0 for c in tb] for i in range(n)]
        d = [sum(a * c[2] for a, c in zip(A[i], tb)) for i in range(n)]
        g = [sum(A[i][j] / d[i] for i in range(n)) for j in range(len(tb))]
        assert max(abs(gj - n) for gj, c in zip(g, tb) if c[2] > 0) / n < 1e-12
        assert all((gj - n) / n <= 1e-12 for gj, c in zip(g, tb) if c[2] == 0)
        assert abs(sum(c[2] for c in tb) - 1) < 1e-14


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


def test_sharing_mfunction_is_proper_and_A_n_on_exact_data():
    x = [1.3, 2.7, 4.1, 6.6, 9.2]
    M = p.npi_mfun(x, x)
    assert len(M) == 6 and all(close(m, 1 / 6, 1e-12) for _, _, m in M)
    assert all(close(a, b, 1e-6) for (a, _, _), b in zip(M, [0.0] + x))
    rng = random.Random(5)
    for _ in range(200):
        n = rng.randint(2, 15)
        a = [round(rng.uniform(1, 10), 2) for _ in range(n)]
        ty = rng.choices("eirl", weights=[0.35, 0.4, 0.15, 0.1], k=n)
        L = [0.0 if t == "l" else x for x, t in zip(a, ty)]
        R = [x if t in "el" else (INF if t == "r" else round(x + rng.uniform(0.3, 7), 2)) for x, t in zip(a, ty)]
        if len(set(a + [r for r, t in zip(R, ty) if t == "i"])) < n + sum(t == "i" for t in ty):
            continue
        M = p.npi_mfun(L, R)
        assert close(sum(m for _, _, m in M), 1.0, 1e-12) and all(lo < hi for lo, hi, _ in M) and len(M) <= 2 * n + 1


def test_nested_regions_are_the_ordering_rule():
    """On nested or disjoint regions the sharing M-function equals the M-function of the ordering rule, in which every
    ordering consistent with the regions is equally likely; the rule is computed here by enumeration."""
    from itertools import permutations

    def rule(L, R):
        n = len(L); acc = {}; E = 0
        for q in permutations(range(n)):
            if any(max(L[i] for i in q[:k]) >= min(R[i] for i in q[k:]) for k in range(1, n)):
                continue
            E += 1
            for k in range(n + 1):
                key = (max((L[i] for i in q[:k]), default=0.0), min((R[i] for i in q[k:]), default=INF))
                acc[key] = acc.get(key, 0) + 1
        return [(a, b, acc[(a, b)] / ((n + 1) * E)) for (a, b) in sorted(acc)]

    def nested(L, R):
        n = len(L)
        for a in range(n):
            for b in range(a + 1, n):
                if L[a] < R[a] and L[b] < R[b] and not ((L[a] <= L[b] and R[b] <= R[a]) or (L[b] <= L[a] and R[a] <= R[b])
                                                        or R[a] < L[b] or R[b] < L[a]):
                    return False
        return True

    rng = random.Random(8); done = 0
    while done < 60:
        n = rng.randint(2, 6)
        a = [round(rng.uniform(1, 10), 2) for _ in range(n)]
        ty = rng.choices("eir", weights=[0.4, 0.45, 0.15], k=n)
        R = [x if t == "e" else (INF if t == "r" else round(x + rng.uniform(0.3, 7), 2)) for x, t in zip(a, ty)]
        ends = a + [r for r, t in zip(R, ty) if t == "i"]
        if len(set(ends)) < len(ends) or not nested(a, R):
            continue
        done += 1
        A = p.npi_mfun(a, R); B = rule(a, R)
        assert len(A) == len(B)
        for (lo, hi, m), (lo2, hi2, m2) in zip(A, B):
            assert abs(lo - lo2) < 1e-5 and (hi == hi2 or abs(hi - hi2) < 1e-5) and close(m, m2, 1e-12)


def test_turnbull_with_right_censoring_far_out():
    """A finite right endpoint far beyond the left endpoints, and evaluation beyond all the data."""
    tb = p.npi_turnbull([0.2, 0.5, 1], [0.2, 20, INF])
    assert [(a, b) for a, b, _ in tb] == [(0.2, 0.2), (1.0, 20.0)]
    assert close(tb[0][2], 1 / 3, 1e-12) and close(tb[1][2], 2 / 3, 1e-12)
    s = p.npi_survival([0.2, 0.5, 1], [0.2, 20, INF], [10, 30], members=True)
    assert close(s[0]["TB.upper"], 2 / 3, 1e-12) and abs(s[1]["TB.upper"]) < 1e-12
    tb = p.npi_turnbull([2, 3], [2, INF])
    assert tb[1][1] == INF and close(tb[1][2], 0.5, 1e-12)
    for r in p.npi_survival([2, 3], [2, INF], [2.5, 100, 1e6], members=True):
        assert close(r["TB.upper"], 0.5, 1e-12) and close(r["P.upper"], 2 / 3, 1e-12) and close(r["upper"], 2 / 3, 1e-12)


def test_ties_one_convention_for_every_member():
    """Half-open reading of a censored region: an event precedes a censoring tied with it, an exact value at the end
    of an interval lies inside it and at the start outside, and a region ending at a value lies below one starting there."""
    L = [1, 1, 1, 1, 1, 2]; R = [1, 1, 1, INF, INF, 2]
    tb = p.npi_turnbull(L, R)
    assert [a for a, _, _ in tb] == [1.0, 2.0] and close(tb[0][2], 0.5, 1e-12) and close(tb[1][2], 0.5, 1e-12)
    s = p.npi_survival(L, R, 1.5, members=True)[0]
    assert close(s["lower"], 2 / 7, 1e-12) and close(s["S.lower"], 2 / 7, 1e-12) and close(s["TB.upper"], 0.5, 1e-12)
    old = p.npi_turnbull(L, R, untie=False)
    assert close(old[0][2], 0.75, 1e-10) and close(old[1][2], 0.25, 1e-10)
    tb = p.npi_turnbull([2, 5], [5, 5])
    assert len(tb) == 1 and tb[0][0] == 5.0 and close(tb[0][2], 1.0, 1e-12)
    tb = p.npi_turnbull([2, 2], [2, 5])
    assert [(a, b) for a, b, _ in tb] == [(2.0, 2.0), (2.0, 5.0)] and close(tb[0][2], 0.5, 1e-12)
    tb = p.npi_turnbull([1, 3, 3], [3, 5, 3])
    assert [(a, b) for a, b, _ in tb] == [(3.0, 3.0), (3.0, 5.0)] and close(tb[0][2], 2 / 3, 1e-12)
    tb = p.npi_turnbull([1, 3], [3, 5])
    assert [(a, b) for a, b, _ in tb] == [(1.0, 3.0), (3.0, 5.0)] and close(tb[0][2], 0.5, 1e-12)


def test_unit_of_measurement_and_origin():
    rng = random.Random(22)
    for _ in range(100):
        n = rng.randint(2, 12)
        a = [rng.randint(1, 8) for _ in range(n)]
        ty = rng.choices("eirl", weights=[0.35, 0.4, 0.15, 0.1], k=n)
        L = [0.0 if t == "l" else float(x) for x, t in zip(a, ty)]
        R = [float(x) if t in "el" else (INF if t == "r" else float(x + rng.randint(1, 4))) for x, t in zip(a, ty)]
        ts = [0.5 + j for j in range(14)]
        s1 = p.npi_survival(L, R, ts, members=True)
        for k in (1e-6, 1e6):
            s2 = p.npi_survival([x * k for x in L], [x * k for x in R], [t * k for t in ts], members=True)
            assert all(abs(r1[key] - r2[key]) < 1e-12 for r1, r2 in zip(s1, s2) for key in r1 if key != "t")
        M = p.npi_mfun(L, R)
        ends = set([0.0] + L + R + [INF])
        assert all(lo in ends and hi in ends for lo, hi, _ in M) and close(sum(m for _, _, m in M), 1.0, 1e-12)
    for bad in ([-2, -1, 3], [-1.0, 2.0]):
        try:
            p.npi_survival(bad, bad, 0)
            assert False
        except ValueError:
            pass
