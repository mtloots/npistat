#!/bin/sh
# Parity harness for npistat: every exported routine on a fixed grid of inputs, byte-identical
# between the R and Python fronts. Both fronts compile the same C source; the guard below refuses
# to run against stale copies.
set -e
OUT=${TMPDIR:-/tmp}/parity_npistat
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/npistat/src" "$PKG/npistat-py/src/npistat" || exit 1
Rscript -e '
suppressMessages(library(npistat))
L <- c(3.5, 3.4, 4.7, 6.9, 5.6, 2.4, 1, 1, 8); R <- c(6.5, 6.8, 5.9, 6.9, 7.9, 6.4, 1, Inf, 8)
tm <- c(3, 4, 6, 7, 8, 9); ev <- c(1, 0, 1, 0, 1, 1)
s <- npi_survival(L, R, c(0.5, 2.25, 5.75, 7.1, 9.5), members = TRUE)
v <- c(as.matrix(s[, -1]), npi_mfun(L, R)$mass, npi_rc(tm, ev)$mass, npi_turnbull(L, R)$mass,
       unlist(npi_bernoulli(1000, 2, 100, 1:5)[, 2:3]), unlist(npi_bernoulli(20000, 0, 100, 1)[, 2:3]),
       unlist(npi_future(L, R, c(2, 6), 10, 2, "above")[, 2:3]), unlist(npi_future(L, R, c(2, 6), 10, 2, "below")[, 2:3]),
       unlist(npi_twolevel(c(0, 0, 0.0035, 0.0068, 0.0061, 0.0017), m = 100, r = 1:4)[, 2:3]),
       { iv <- npi_rates(c(0, 2, 1, 0, 3), c(1031, 1062, 1083, 1500, 2900)); unlist(npi_twolevel(iv$L, iv$R, 3539, 1:3)[, 2:3]) })
cat(sprintf("%.17g", v), sep = "\n")
' > "$OUT/r.txt"
PYTHONPATH="$PKG/npistat-py/src" python3 - <<PY > "$OUT/py.txt"
import npistat as p
INF = float("inf")
L = [3.5, 3.4, 4.7, 6.9, 5.6, 2.4, 1, 1, 8]; R = [6.5, 6.8, 5.9, 6.9, 7.9, 6.4, 1, INF, 8]
tm = [3, 4, 6, 7, 8, 9]; ev = [1, 0, 1, 0, 1, 1]
s = p.npi_survival(L, R, [0.5, 2.25, 5.75, 7.1, 9.5], members=True)
keys = ["lower", "upper", "S.lower", "S.upper", "P.lower", "P.upper", "C.lower", "C.upper", "TB.lower", "TB.upper"]
v = [r[k] for k in keys for r in s]
v += [m[2] for m in p.npi_mfun(L, R)] + [m[2] for m in p.npi_rc(tm, ev)] + [c[2] for c in p.npi_turnbull(L, R)]
b = p.npi_bernoulli(1000, 2, 100, [1, 2, 3, 4, 5]); v += [x[1] for x in b] + [x[2] for x in b]
b = p.npi_bernoulli(20000, 0, 100, 1); v += [b[0][1], b[0][2]]
for side in ("above", "below"):
    f = p.npi_future(L, R, [2, 6], 10, 2, side); v += [x[1] for x in f] + [x[2] for x in f]
t2 = p.npi_twolevel([0, 0, 0.0035, 0.0068, 0.0061, 0.0017], None, 100, [1, 2, 3, 4]); v += [x[1] for x in t2] + [x[2] for x in t2]
a, b = p.npi_rates([0, 2, 1, 0, 3], [1031, 1062, 1083, 1500, 2900])
t2 = p.npi_twolevel(a, b, 3539, [1, 2, 3]); v += [x[1] for x in t2] + [x[2] for x in t2]
print("\n".join("%.17g" % x for x in v))
PY
if cmp -s "$OUT/r.txt" "$OUT/py.txt"; then
  echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "VERDICT: MISMATCH"; diff "$OUT/r.txt" "$OUT/py.txt" | head; exit 1
fi
