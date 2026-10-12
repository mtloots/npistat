"""npistat: Python (ctypes) binding to the shared C back-end. Pure standard library.

The same C sources back the R package 'npistat'; the parity harness in the repository requires every
exported quantity to agree between the two fronts, which holds because every computation is in C.

LOADING has two paths, as in the other packages in this programme. An INSTALLED package carries a
library that pip built at install time, so no compiler is needed here. A SOURCE CHECKOUT has the
.c file but no library, so it is compiled on demand into a writable cache directory.
"""
import ctypes as _ct
import glob as _glob
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))


def _find_built():
    for pat in ("_libnpistat*.so", "_libnpistat*.dylib", "_libnpistat*.pyd",
                "_libnpistat*.dll"):
        hits = sorted(_glob.glob(_os.path.join(_here, pat)))
        if hits:
            return hits[0]
    return None


def _cache_dir():
    if _os.access(_here, _os.W_OK):
        return _here
    base = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(_os.path.expanduser("~"), ".cache")
    d = _os.path.join(base, "npistat")
    _os.makedirs(d, exist_ok=True)
    return d


def _ext():
    return "dll" if _sys.platform == "win32" else ("dylib" if _sys.platform == "darwin" else "so")


def _compile_from_source():
    import subprocess as _sp
    srcs = sorted(f for f in _glob.glob(_os.path.join(_here, "*.c"))
                  if _os.path.basename(f) != "init.c")
    if not srcs:
        return None
    out = _os.path.join(_cache_dir(), "libnpistat." + _ext())
    cc = _os.environ.get("CC") or "cc"
    try:
        _sp.check_call([cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out])
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "npistat could not build its C back-end from the bundled sources (%s). In an "
            "installed copy this should never happen, because pip builds the library at install "
            "time. In a source checkout a C compiler is required; set CC if it is not 'cc'." % exc)
    return out



_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libnpistat." + _ext())
    _libpath = _cached if _os.path.exists(_cached) else _compile_from_source()
if _libpath is None:
    raise ImportError("npistat: no compiled back-end and no C sources to build one from.")

_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)

for _f, _sig in (
    ("C_npi_bern", [_ip, _ip, _ip, _ip, _ip, _dp, _dp]),
    ("C_npi_mfuture", [_ip, _dp, _dp, _ip, _ip, _ip, _ip, _dp, _dp, _dp]),
    ("C_npi_untie", [_ip, _dp, _dp, _ip, _dp, _ip, _dp, _dp, _dp]),
    ("C_npi_share", [_ip, _dp, _dp, _ip, _dp]),
    ("C_npi_place", [_ip, _dp, _dp, _ip, _dp, _dp, _dp, _dp]),
    ("C_npi_mbounds", [_ip, _dp, _dp, _dp, _ip, _ip, _dp, _dp, _dp]),
    ("C_npi_turnbull", [_ip, _dp, _dp, _ip, _dp, _ip, _dp, _dp, _dp]),
    ("C_npi_bounds", [_ip, _dp, _dp, _ip, _dp, _ip, _ip, _dp] + [_dp] * 10),
    ("C_npi_rc", [_ip, _dp, _ip, _dp, _dp, _dp]),
    ("C_npi_twolevel", [_ip, _dp, _dp, _ip, _ip, _ip, _ip, _dp, _dp]),
):
    getattr(_lib, _f).argtypes = _sig
    getattr(_lib, _f).restype = None

_INF = float("inf")


def _d(xs):
    xs = list(xs)
    return (_ct.c_double * max(len(xs), 1))(*[float(x) for x in xs])


def _i(xs):
    xs = list(xs)
    return (_ct.c_int * max(len(xs), 1))(*[int(x) for x in xs])


def _z(k):
    return (_ct.c_double * max(int(k), 1))()


def _one(x):
    return _ct.byref(_ct.c_int(int(x)))


def _oned(x):
    return _ct.byref(_ct.c_double(float(x)))


def _vec(x):
    return [float(v) for v in (x if hasattr(x, "__len__") else [x])]


def _regions(L, R):
    L, R = _vec(L), _vec(R)
    if len(L) != len(R):
        raise ValueError("L and R must have the same length")
    if any(v != v for v in L + R):
        raise ValueError("L and R must not contain NaN")
    if any(r < l for l, r in zip(L, R)):
        raise ValueError("every region needs L <= R")
    if any(abs(l) == _INF for l in L):
        raise ValueError("L must be finite; code a left-censored observation by L = 0")
    return L, R, [int(l == r) for l, r in zip(L, R)], len(L)


def npi_bernoulli(n, k, m, r=1):
    """Lower and upper probabilities that at least r of m future observations are events, given k
    events among n exchangeable past observations (Coolen 1998). r may be a list.

    Returns a list of (r, lower, upper)."""
    n, k, m = int(n), int(k), int(m)
    if n < 1 or k < 0 or k > n or m < 1:
        raise ValueError("need n >= 1, 0 <= k <= n and m >= 1")
    rs = [int(v) for v in (r if hasattr(r, "__len__") else [r])]
    lo, up = _z(len(rs)), _z(len(rs))
    _lib.C_npi_bern(_one(n), _one(k), _one(m), _i(rs), _one(len(rs)), lo, up)
    return [(rs[j], lo[j], up[j]) for j in range(len(rs))]


def npi_future(L, R, t, m, r=1, side="above"):
    """Bounds that at least r of m future observations fall above t (side="above") or at or below
    it (side="below"), from regions [L_i, R_i] by completions. Returns a list of (t, lower, upper)."""
    if side not in ("above", "below"):
        raise ValueError("side must be 'above' or 'below'")
    L, R, ex, n = _regions(L, R)
    ts = _vec(t)
    lo, up = _z(len(ts)), _z(len(ts))
    _lib.C_npi_mfuture(_one(n), _d(L), _d(R), _one(m), _one(r), _one(side == "below"),
                       _one(len(ts)), _d(ts), lo, up)
    return [(ts[j], lo[j], up[j]) for j in range(len(ts))]


def npi_mfun(L, R, closedR=True, eps=1e-7):
    """M-functions of the sharing specification M_S for regions [L_i, R_i] (R_i = inf for right
    censoring, L_i == R_i for an exact value). Returns a list of (from, to, mass), merged where an
    interval occurs twice; under right censoring these are the rc-A(n) values of Coolen and Yan."""
    L, R, ex, n = _regions(L, R)
    L2, R2 = _z(n), _z(n)
    _lib.C_npi_untie(_one(n), _d(L), _d(R), _i(ex), _oned(eps), _one(bool(closedR)), _oned(1e-9), L2, R2)
    own = _z(n + 1)
    _lib.C_npi_share(_one(n), L2, R2, _i(ex), own)
    lo, hi, mm = _z(n + 1), _z(n + 1), _z(n + 1)
    _lib.C_npi_place(_one(n), L2, R2, _i(ex), own, lo, hi, mm)
    acc = {}
    for j in range(n + 1):
        key = (round(lo[j], 6), round(hi[j], 6) if hi[j] != _INF else _INF)
        acc[key] = acc.get(key, 0.0) + mm[j]
    return [(a, b, acc[(a, b)]) for (a, b) in sorted(acc)]


def npi_rc(time, status):
    """The rc-A(n) M-functions of Coolen and Yan (2004) in closed form, no ties; status 1 is an
    event, 0 a right censoring. Returns a list of (from, to, mass) in time order."""
    time, status = _vec(time), [int(s) for s in status]
    n = len(time)
    if len(status) != n:
        raise ValueError("time and status must have the same length")
    if len(set(time)) != n:
        raise ValueError("the closed form assumes no ties; use npi_mfun()")
    lo, hi, mm = _z(n + 1), _z(n + 1), _z(n + 1)
    _lib.C_npi_rc(_one(n), _d(time), _i(status), lo, hi, mm)
    return [(lo[j], hi[j], mm[j]) for j in range(n + 1)]


def npi_turnbull(L, R, maxit=20000, tol=1e-13):
    """The Turnbull estimator: maximal cliques and masses by self-consistency. Returns a list of
    (left, right, mass); an infinite right endpoint is reported as 3 max(L) + 1."""
    L, R, ex, n = _regions(L, R)
    K = _ct.c_int(0)
    ql, qr, p = _z(n), _z(n), _z(n)
    _lib.C_npi_turnbull(_one(n), _d(L), _d(R), _one(maxit), _oned(tol), _ct.byref(K), ql, qr, p)
    return [(ql[j], qr[j], p[j]) for j in range(K.value)]


def npi_survival(L, R, t, closedR=True, eps=1e-7, members=False):
    """Lower and upper survival probabilities P(X_{n+1} > t) of the construction: the union of the
    sharing and anchored specifications, clipped to the completion bound. Evaluate away from the
    data values. Returns a list of dicts with t, lower, upper and, if members, the bounds of
    M_S, M_P, M_C and the Turnbull set."""
    L, R, ex, n = _regions(L, R)
    ts = _vec(t)
    nt = len(ts)
    out = [_z(nt) for _ in range(10)]
    _lib.C_npi_bounds(_one(n), _d(L), _d(R), _i(ex), _oned(eps), _one(bool(closedR)), _one(nt),
                      _d(ts), *out)
    names = ["lower", "upper", "S.lower", "S.upper", "P.lower", "P.upper", "C.lower", "C.upper",
             "TB.lower", "TB.upper"]
    res = []
    for j in range(nt):
        row = {"t": ts[j]}
        for k, nm in enumerate(names if members else names[:2]):
            row[nm] = out[k][j]
        res.append(row)
    return res


def npi_twolevel(L, R=None, m=1, r=1, closedR=True):
    """Two-level bounds: lower and upper probabilities that at least r of m units have an event in
    the next period, when past periods are exchangeable (Hill's A(n) on their rates, given as
    proportions; L == R for a rate taken as known, L < R for an interval) and units within a period
    are exchangeable (binomial given the period's rate). The upper probability is always at least
    1/(n+1). Returns a list of (r, lower, upper)."""
    L = _vec(L)
    R = L[:] if R is None else _vec(R)
    if len(L) != len(R):
        raise ValueError("L and R must have the same length")
    if any(not (0.0 <= l <= u <= 1.0) for l, u in zip(L, R)):
        raise ValueError("rates must satisfy 0 <= L <= R <= 1")
    rs = [int(v) for v in (r if hasattr(r, "__len__") else [r])]
    lo, up = _z(len(rs)), _z(len(rs))
    _lib.C_npi_twolevel(_one(len(L)), _d(L), _d(R), _one(m), _i(rs), _one(len(rs)), _one(bool(closedR)), lo, up)
    return [(rs[j], lo[j], up[j]) for j in range(len(rs))]


def npi_rates(k, N):
    """Interval-valued period rates [k/(N+1), (k+1)/(N+1)] from k events among N units in each
    period. Returns (L, R) as two lists."""
    k, N = _vec(k), _vec(N)
    if len(k) != len(N) or any(a < 0 or a > b for a, b in zip(k, N)):
        raise ValueError("need 0 <= k <= N, elementwise")
    return [a / (b + 1) for a, b in zip(k, N)], [(a + 1) / (b + 1) for a, b in zip(k, N)]
