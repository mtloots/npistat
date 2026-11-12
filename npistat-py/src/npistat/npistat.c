/* npistat: nonparametric predictive inference (NPI) for censored and Bernoulli data.
 *
 * One pure-C back-end, libm only, shared byte for byte by the R package (.C interface, init.c)
 * and the Python package (ctypes). Every argument is a pointer, as the .C ABI requires.
 *
 * Two families of computation.
 *
 * 1. Censored data (Loots, "Nonparametric predictive inference for interval-censored data:
 *    M-functions anchored to the self-consistent estimator"). Each observation is known only to
 *    lie in a region with endpoints L_i <= R_i; R_i = +infinity for right censoring, L_i = R_i for
 *    an exact observation. Ties are broken first (C_npi_untie), a censored region being read as
 *    half-open, (L_i, R_i], and every specification is computed from the regions so obtained:
 *      M_S  the sharing assumption ic-A(n), which is Coolen and Yan's rc-A(n) under right censoring;
 *      M_P  n/(n+1) of each Turnbull clique mass on its clique, 1/(n+1) on the whole support;
 *      M_C  the completion bound, the envelope of A(n) over all completions of the data.
 *    The construction's survival bounds are min(S_lo, P_lo) and max(S_up, P_up); both members lie
 *    inside M_C.
 *
 * 2. Bernoulli data and m future observations (Coolen 1998; Coolen and Coolen-Schrijner 2005),
 *    by counting allocations of m future observations to the n+1 intervals of A(n), and the same
 *    count applied to censored data through completions.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------------------------- */
/* allocation counting                                                                          */
/* ------------------------------------------------------------------------------------------- */

static double lchoose_(double a, double b) { return lgamma(a + 1.0) - lgamma(b + 1.0) - lgamma(a - b + 1.0); }

/* log of the number of ways to put k indistinguishable future observations into `boxes`
   intervals; with no intervals the count is 1 when k == 0 and 0 otherwise */
static double lbars(int k, int boxes) {
    if (boxes == 0) return k == 0 ? 0.0 : -INFINITY;
    return lchoose_((double)(k + boxes - 1), (double)(boxes - 1));
}

/* Lower probability that at least r of m future observations fall strictly above a threshold
   when i of the n past observations lie at or below it (Coolen 1998). */
static double plo_(int m, int n, int i, int r) {
    if (r > m) return 0.0;
    if (r <= 0) return 1.0;
    if (i < 0) return 1.0;
    if (i > n) return 0.0;
    double lt = lchoose_((double)(n + m), (double)m), s = 0.0;
    for (int d = r; d <= m; d++) {
        double a = lbars(d, n - i), b = lbars(m - d, i + 1);
        if (isinf(a) || isinf(b)) continue;
        s += exp(a + b - lt);
    }
    return s > 1.0 ? 1.0 : s;
}
/* the upper probability is the lower one with the straddling interval counted above */
static double pup_(int m, int n, int i, int r) { return plo_(m, n, i - 1, r); }

/* 1 - plo_(m, n, i, r), summed directly over d = 0, ..., r - 1 rather than by complement: the
   complement of a probability near one loses relative accuracy in proportion to 1/(1 - p), which
   at n = 20000 and m = 100 costs eight digits. The terms of the counting identity sum to one, so
   the direct sum is the same quantity without the cancellation. */
static double qlo_(int m, int n, int i, int r) {
    if (r > m) return 1.0;
    if (r <= 0) return 0.0;
    if (i < 0) return 0.0;
    if (i > n) return 1.0;
    double lt = lchoose_((double)(n + m), (double)m), s = 0.0;
    for (int d = 0; d < r; d++) {
        double a = lbars(d, n - i), b = lbars(m - d, i + 1);
        if (isinf(a) || isinf(b)) continue;
        s += exp(a + b - lt);
    }
    return s > 1.0 ? 1.0 : s;
}
static double qup_(int m, int n, int i, int r) { return qlo_(m, n, i - 1, r); }

/* Bernoulli data: k events among n, bounds for at least r events among m future, r vectorised.
   An event is an observation below the threshold, so the bounds follow by complementation. */
void C_npi_bern(const int *n, const int *k, const int *m, const int *r, const int *nr,
                double *lo, double *up) {
    for (int j = 0; j < *nr; j++) {
        lo[j] = qup_(*m, *n, *k, *m - r[j] + 1);
        up[j] = qlo_(*m, *n, *k, *m - r[j] + 1);
    }
}

/* Censored data, m future observations, by completions: bounds for at least r of m future
   observations above t (side = 0) or at or below t (side = 1), t vectorised. */
void C_npi_mfuture(const int *n, const double *L, const double *R, const int *m, const int *r,
                   const int *side, const int *nt, const double *t, double *lo, double *up) {
    for (int j = 0; j < *nt; j++) {
        int iL = 0, iR = 0;
        for (int i = 0; i < *n; i++) { if (L[i] <= t[j]) iL++; if (R[i] <= t[j]) iR++; }
        if (*side == 0) {
            lo[j] = plo_(*m, *n, iL, *r);
            up[j] = pup_(*m, *n, iR, *r);
        } else {
            lo[j] = qup_(*m, *n, iR, *m - *r + 1);
            up[j] = qlo_(*m, *n, iL, *m - *r + 1);
        }
    }
}

/* ------------------------------------------------------------------------------------------- */
/* ties                                                                                         */
/* ------------------------------------------------------------------------------------------- */

static const double *sL_, *sR_; static const int *sE_;
static int cmp_untie(const void *a, const void *b) {
    int i = *(const int *)a, j = *(const int *)b;
    if (sL_[i] != sL_[j]) return sL_[i] < sL_[j] ? -1 : 1;
    int ci = !sE_[i], cj = !sE_[j];                      /* events before censorings */
    if (ci != cj) return ci < cj ? -1 : 1;
    if (sR_[i] != sR_[j]) return sR_[i] > sR_[j] ? -1 : 1; /* containing region first */
    return i < j ? -1 : (i > j);
}
static int cmp_dbl(const void *a, const void *b);
/* the largest of the m increasing values v that is at most x */
static double floor_of(const double *v, int m, double x) {
    int lo = 0, hi = m - 1;
    while (lo < hi) { int mid = (lo + hi + 1) / 2; if (v[mid] <= x) lo = mid; else hi = mid - 1; }
    return v[lo];
}

/* Step 1, the tie convention: Hill's H(n) by the epsilon limit (Coolen and Yan 2004, Section 7). A
   censored region is read as half-open, (L_i, R_i], an exact observation as the point L_i = R_i, and
   ties are broken by increments that stand for infinitesimals:
     - endpoint values closer than tol times the largest finite endpoint are taken as tied;
     - the left endpoints are taken in increasing order, exact observations before censored ones at
       a tied value and the containing region first among tied censored ones, and the k-th of them
       receives the increment k eps;
     - when closedR is 1 the right endpoint of a censored region is moved just past the exact
       observations tied with it, and no further, so that those observations lie inside the region
       and the region still ends before any censored region that starts at the same value; when
       closedR is 0 right endpoints are not moved, and an exact observation tied with one lies
       outside. Tied right endpoints stay tied.
   eps <= 0 asks for an increment chosen from the data, the smallest gap between distinct endpoint
   values divided by 64 (n + 2): no endpoint then moves by more than a sixty-fourth of that gap, so
   the order of distinct values is never changed and the result does not depend on the unit of
   measurement. */
void C_npi_untie(const int *n, const double *L, const double *R, const int *ex, const double *eps,
                 const int *closedR, const double *tol, double *L2, double *R2) {
    int nn = *n, m = 0, mc = 0;
    if (nn <= 0) return;
    int *o = (int *)malloc(sizeof(int) * nn);
    double *Lc = (double *)malloc(sizeof(double) * nn), *Rc = (double *)malloc(sizeof(double) * nn);
    double *v = (double *)malloc(sizeof(double) * 2 * nn);
    double scale = 0.0, gap = INFINITY;
    for (int i = 0; i < nn; i++) { v[m++] = L[i]; if (isfinite(R[i])) v[m++] = R[i]; }
    for (int j = 0; j < m; j++) if (fabs(v[j]) > scale) scale = fabs(v[j]);
    if (!(scale > 0.0)) scale = 1.0;
    qsort(v, m, sizeof(double), cmp_dbl);
    double thr = (*tol > 0.0 ? *tol : 0.0) * scale;
    for (int j = 0; j < m; j++)                           /* one representative for each set of tied values */
        if (mc == 0 || v[j] - v[mc - 1] > thr) { if (mc > 0 && v[j] - v[mc - 1] < gap) gap = v[j] - v[mc - 1]; v[mc++] = v[j]; }
    if (!isfinite(gap)) gap = scale;
    double e = *eps > 0.0 ? *eps : gap / (64.0 * (nn + 2));
    for (int i = 0; i < nn; i++) { Lc[i] = floor_of(v, mc, L[i]); Rc[i] = isfinite(R[i]) ? floor_of(v, mc, R[i]) : R[i]; o[i] = i; }
    sL_ = Lc; sR_ = Rc; sE_ = ex;
    qsort(o, nn, sizeof(int), cmp_untie);
    for (int k = 0; k < nn; k++) L2[o[k]] = Lc[o[k]] + e * (double)(k + 1);
    for (int i = 0; i < nn; i++) {
        if (ex[i] || Rc[i] == Lc[i]) { R2[i] = L2[i]; continue; }
        R2[i] = Rc[i];
        if (*closedR && isfinite(Rc[i])) {
            double top = -INFINITY;
            for (int j = 0; j < nn; j++) if (ex[j] && Lc[j] == Rc[i] && L2[j] > top) top = L2[j];
            if (isfinite(top)) R2[i] = top + 0.5 * e;
        }
    }
    free(o); free(Lc); free(Rc); free(v);
}

/* ------------------------------------------------------------------------------------------- */
/* the sharing specification M_S                                                                */
/* ------------------------------------------------------------------------------------------- */

/* compatible set C_i = { j : L_j >= L_i and R_j <= R_i }, itself included */
static int in_C(const double *L, const double *R, int i, int j) { return L[j] >= L[i] && R[j] <= R[i]; }

/* ic-A(n): censored observations share their mass equally over C_i, in decreasing order of |C_i|,
   those of equal size from a common snapshot. own has n+1 entries: own[0] is the region below
   the smallest observation, own[i+1] the mass of observation i. */
void C_npi_share(const int *n, const double *L, const double *R, const int *ex, double *own) {
    int nn = *n;
    int *sz = (int *)calloc(nn > 0 ? nn : 1, sizeof(int));
    int *done = (int *)calloc(nn > 0 ? nn : 1, sizeof(int));
    double *snap = (double *)malloc(sizeof(double) * (nn + 1));
    double *d = (double *)malloc(sizeof(double) * (nn + 1));
    for (int i = 0; i <= nn; i++) own[i] = 1.0 / (nn + 1);
    for (int i = 0; i < nn; i++) {
        if (ex[i]) { done[i] = 1; continue; }
        for (int j = 0; j < nn; j++) if (in_C(L, R, i, j)) sz[i]++;
    }
    for (;;) {
        int s = -1;
        for (int i = 0; i < nn; i++) if (!done[i] && sz[i] > s) s = sz[i];
        if (s < 0) break;
        memcpy(snap, own, sizeof(double) * (nn + 1));
        memset(d, 0, sizeof(double) * (nn + 1));
        for (int i = 0; i < nn; i++) {
            if (done[i] || sz[i] != s) continue;
            double sh = snap[i + 1] / (double)s;
            d[i + 1] += sh - snap[i + 1];
            for (int j = 0; j < nn; j++) if (j != i && in_C(L, R, i, j)) d[j + 1] += sh;
        }
        for (int i = 0; i < nn; i++) if (!done[i] && sz[i] == s) done[i] = 1;
        for (int i = 0; i <= nn; i++) own[i] += d[i];
    }
    free(sz); free(done); free(snap); free(d);
}

/* Place the masses on intervals. Each mass has a left end and a right end: observation i starts at
   L_i and ends at R_i, the unit below the smallest observation starts at the origin, and one unit
   ends at infinity. The left ends and the right ends are paired in increasing order, mass for mass,
   which gives at most 2n+1 intervals; K receives their number, and lo, hi, m need 2n+1 entries.
   On exact data these are the intervals of A(n), on right-censored data those of rc-A(n). */
static const double *pv_;
static int cmp_place(const void *a, const void *b) {
    int i = *(const int *)a, j = *(const int *)b;
    if (pv_[i] != pv_[j]) return pv_[i] < pv_[j] ? -1 : 1;
    return i < j ? -1 : (i > j);
}
void C_npi_place(const int *n, const double *L, const double *R, const int *ex, const double *own,
                 int *K, double *lo, double *hi, double *m) {
    int nn = *n, k = 0, a = 0, b = 0;
    const double tiny = 1e-14;
    double *va = (double *)malloc(sizeof(double) * (nn + 1)), *vb = (double *)malloc(sizeof(double) * (nn + 1));
    double *ma = (double *)malloc(sizeof(double) * (nn + 1)), *mb = (double *)malloc(sizeof(double) * (nn + 1));
    int *oa = (int *)malloc(sizeof(int) * (nn + 1)), *ob = (int *)malloc(sizeof(int) * (nn + 1));
    (void)ex;
    va[0] = -INFINITY; ma[0] = own[0];                    /* the origin, first among the left ends */
    for (int i = 0; i < nn; i++) { va[i + 1] = L[i]; ma[i + 1] = own[i + 1]; }
    for (int i = 0; i < nn; i++) { vb[i] = R[i]; mb[i] = own[i + 1]; }
    vb[nn] = INFINITY; mb[nn] = own[0];                   /* the unit that ends at infinity */
    for (int i = 0; i <= nn; i++) { oa[i] = i; ob[i] = i; }
    pv_ = va; qsort(oa, nn + 1, sizeof(int), cmp_place);
    pv_ = vb; qsort(ob, nn + 1, sizeof(int), cmp_place);
    double ra = ma[oa[0]], rb = mb[ob[0]];
    while (a <= nn && b <= nn) {
        double piece = ra < rb ? ra : rb;
        if (piece > tiny) { lo[k] = oa[a] == 0 ? 0.0 : va[oa[a]]; hi[k] = vb[ob[b]]; m[k] = piece; k++; }
        ra -= piece; rb -= piece;
        if (ra <= tiny) { a++; if (a <= nn) ra += ma[oa[a]]; }
        if (rb <= tiny) { b++; if (b <= nn) rb += mb[ob[b]]; }
    }
    *K = k;
    free(va); free(vb); free(ma); free(mb); free(oa); free(ob);
}

/* Lower and upper probability of (t, infinity) for an M-function on K intervals: the mass on
   intervals inside the event and the mass on intervals meeting it. closed = 0 for open intervals
   (lower counts lo >= t), closed = 1 for closed cliques (lower counts lo > t). */
void C_npi_mbounds(const int *K, const double *lo, const double *hi, const double *m,
                   const int *closed, const int *nt, const double *t, double *low, double *upp) {
    for (int j = 0; j < *nt; j++) {
        double a = 0.0, b = 0.0;
        for (int k = 0; k < *K; k++) {
            if (*closed ? lo[k] > t[j] : lo[k] >= t[j]) a += m[k];
            if (hi[k] > t[j]) b += m[k];
        }
        low[j] = a; upp[j] = b;
    }
}

/* ------------------------------------------------------------------------------------------- */
/* the Turnbull estimator                                                                       */
/* ------------------------------------------------------------------------------------------- */

static int cmp_dbl(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y);
}
static int uniq(double *v, int k) {
    qsort(v, k, sizeof(double), cmp_dbl);
    int u = 0;
    for (int i = 0; i < k; i++) if (u == 0 || v[i] != v[u - 1]) v[u++] = v[i];
    return u;
}

/* The NPMLE by a damped active-set Newton method. On a support S the masses x maximise the concave
   phi(x) = sum_i log d_i - n sum_{j in S} x_j, d_i = sum_{j in S} a_ij x_j, whose maximiser has sum x = 1 and
   satisfies the optimality conditions sum_i a_ij / d_i = n on S. Each Newton step is limited to keep the masses
   positive: a mass that reaches zero leaves S, and a backtracking search keeps phi increasing. Once settled,
   the conditions sum_i a_ij / d_i <= n are checked off S and the worst violator joins S. Returns 1 when the
   conditions hold to machine precision (p then holds the masses), 0 on a singular system or no convergence. */
static int tb_inC(const double *L, const double *R, const double *ql, const double *qr, int i, int j) { return ql[j] >= L[i] && qr[j] <= R[i]; }
static double tb_phi(int nn, int ns, const int *S, const double *x, const double *L, const double *R, const double *ql, const double *qr) {
    double f = 0.0, sx = 0.0;
    for (int c = 0; c < ns; c++) sx += x[c];
    for (int i = 0; i < nn; i++) { double s = 0.0; for (int c = 0; c < ns; c++) if (tb_inC(L, R, ql, qr, i, S[c])) s += x[c]; if (!(s > 0.0)) return -INFINITY; f += log(s); }
    return f - nn * sx;
}
static int tb_newton(int nn, int k, const double *L, const double *R, const double *ql, const double *qr, double *p) {
    int *S = (int *)malloc(sizeof(int) * k), *inS = (int *)calloc(k, sizeof(int)), ns = 0, ok = 0;
    double *x = (double *)malloc(sizeof(double) * k), *xt = (double *)malloc(sizeof(double) * k), *d = (double *)malloc(sizeof(double) * nn);
    double *F = (double *)malloc(sizeof(double) * k), *J = (double *)malloc(sizeof(double) * k * k), *D = (double *)malloc(sizeof(double) * k);
    double pm = 0.0;
    for (int j = 0; j < k; j++) if (p[j] > pm) pm = p[j];
    for (int j = 0; j < k; j++) if (p[j] > 1e-10 * pm) { S[ns] = j; x[ns] = p[j]; inS[j] = 1; ns++; }
    for (int outer = 0; outer < 4 * k + 16 && !ok; outer++) {
        int settled = 0;
        for (int it = 0; it < 400; it++) {
            double xm0 = 0.0; for (int c = 0; c < ns; c++) if (x[c] > xm0) xm0 = x[c];
            for (int c = 0; c < ns && ns > 1; c++) if (x[c] < 1e-14 * xm0) {     /* a negligible mass leaves the support */
                inS[S[c]] = 0; for (int q = c; q + 1 < ns; q++) { S[q] = S[q + 1]; x[q] = x[q + 1]; } ns--; c--; }
            for (int i = 0; i < nn; i++) { double s = 0.0; for (int c = 0; c < ns; c++) if (tb_inC(L, R, ql, qr, i, S[c])) s += x[c]; d[i] = s; if (!(s > 0.0)) goto done; }
            for (int a = 0; a < ns; a++) {
                double f = -(double)nn;
                for (int i = 0; i < nn; i++) if (tb_inC(L, R, ql, qr, i, S[a])) f += 1.0 / d[i];
                F[a] = f;
                for (int b = a; b < ns; b++) {
                    double s = 0.0;
                    for (int i = 0; i < nn; i++) if (tb_inC(L, R, ql, qr, i, S[a]) && tb_inC(L, R, ql, qr, i, S[b])) s -= 1.0 / (d[i] * d[i]);
                    J[a * ns + b] = s; J[b * ns + a] = s;
                }
            }
            for (int a = 0; a < ns; a++) D[a] = -F[a];                    /* solve J D = -F, partial pivoting */
            double scale = 0.0; for (int q = 0; q < ns * ns; q++) if (fabs(J[q]) > scale) scale = fabs(J[q]);
            for (int c = 0; c < ns; c++) {
                int piv = c; double big = fabs(J[c * ns + c]);
                for (int r = c + 1; r < ns; r++) if (fabs(J[r * ns + c]) > big) { big = fabs(J[r * ns + c]); piv = r; }
                if (!(big > 1e-14 * scale)) goto done;                     /* singular: the masses are not identified */
                if (piv != c) { for (int q = 0; q < ns; q++) { double t = J[c * ns + q]; J[c * ns + q] = J[piv * ns + q]; J[piv * ns + q] = t; } double t = D[c]; D[c] = D[piv]; D[piv] = t; }
                for (int r = c + 1; r < ns; r++) { double f = J[r * ns + c] / J[c * ns + c]; if (f != 0.0) { for (int q = c; q < ns; q++) J[r * ns + q] -= f * J[c * ns + q]; D[r] -= f * D[c]; } }
            }
            for (int c = ns - 1; c >= 0; c--) { double s = D[c]; for (int q = c + 1; q < ns; q++) s -= J[c * ns + q] * D[q]; D[c] = s / J[c * ns + c]; }
            double amax = INFINITY; int cz = -1;
            for (int c = 0; c < ns; c++) if (D[c] < 0.0) { double a = -x[c] / D[c]; if (a < amax) { amax = a; cz = c; } }
            double f0 = tb_phi(nn, ns, S, x, L, R, ql, qr), slope = 0.0;
            for (int c = 0; c < ns; c++) slope += F[c] * D[c];
            if (amax <= 1.0 && ns > 1) {                                  /* a mass reaches zero within the full step */
                for (int c = 0; c < ns; c++) xt[c] = x[c] + amax * D[c];
                xt[cz] = 0.0;
                if (tb_phi(nn, ns, S, xt, L, R, ql, qr) >= f0 - 1e-12 * fabs(f0)) {
                    inS[S[cz]] = 0; for (int c = cz; c + 1 < ns; c++) { S[c] = S[c + 1]; xt[c] = xt[c + 1]; } ns--;
                    for (int c = 0; c < ns; c++) x[c] = xt[c] > 0.0 ? xt[c] : 0.0;
                    continue;
                }
            }
            double fmax = 0.0; for (int c = 0; c < ns; c++) if (fabs(F[c]) > fmax) fmax = fabs(F[c]);
            double a = amax <= 1.0 ? 0.5 * amax : 1.0;
            if (!(amax > 1.0 && fmax < 1e-6 * nn)) {                       /* far from the optimum: backtrack for ascent */
                for (;;) {
                    for (int c = 0; c < ns; c++) xt[c] = x[c] + a * D[c];
                    double f1 = tb_phi(nn, ns, S, xt, L, R, ql, qr);
                    if (f1 >= f0 + 1e-4 * a * slope || a < 1e-12) break;
                    a *= 0.5;
                }
                if (a < 1e-12) break;                                     /* no ascent: leave the decision to the residual */
            } else for (int c = 0; c < ns; c++) xt[c] = x[c] + D[c];      /* near it: full Newton steps, quadratic convergence */
            double step = 0.0, xm = 0.0;
            for (int c = 0; c < ns; c++) { if (fabs(a * D[c]) > step) step = fabs(a * D[c]); x[c] = xt[c]; if (x[c] > xm) xm = x[c]; }
            if (a == 1.0 && (step <= 1e-15 * xm || fmax <= 1e-13 * nn)) { settled = 1; break; }
        }
        for (int i = 0; i < nn; i++) { double s = 0.0; for (int c = 0; c < ns; c++) if (tb_inC(L, R, ql, qr, i, S[c])) s += x[c]; d[i] = s; if (!(s > 0.0)) goto done; }
        double res = 0.0;                                                 /* the equations on the support must hold */
        for (int c = 0; c < ns; c++) { double g = 0.0; for (int i = 0; i < nn; i++) if (tb_inC(L, R, ql, qr, i, S[c])) g += 1.0 / d[i]; if (fabs(g - nn) > res) res = fabs(g - nn); }
        if (!settled && res > 1e-10 * nn) goto done;
        int add = -1; double gv = (double)nn * (1.0 + 1e-12);          /* conditions off the support */
        for (int j = 0; j < k; j++) if (!inS[j]) {
            double g = 0.0; for (int i = 0; i < nn; i++) if (tb_inC(L, R, ql, qr, i, j)) g += 1.0 / d[i];
            if (g > gv) { gv = g; add = j; }
        }
        if (add < 0) { ok = 1; break; }
        double xmin = 1.0; for (int c = 0; c < ns; c++) if (x[c] < xmin) xmin = x[c];
        S[ns] = add; x[ns] = 1e-3 * xmin; inS[add] = 1; ns++;
    }
done:
    if (ok) { double tot = 0.0; for (int c = 0; c < ns; c++) tot += x[c]; for (int j = 0; j < k; j++) p[j] = 0.0; for (int c = 0; c < ns; c++) p[S[c]] = x[c] / tot; }
    free(S); free(inS); free(x); free(xt); free(d); free(F); free(J); free(D);
    return ok;
}

/* Maximal cliques (innermost intervals [l, r] with no other endpoint strictly inside) and the NPMLE of
   Turnbull (1976). The self-consistency iteration from the uniform start locates the support, and
   Newton's method on the optimality conditions (tb_newton) then fixes the masses to machine precision;
   the iteration alone can stop well short of the optimum, notably on current-status data. When the
   masses are not identified (singular Newton system) the iteration runs on to maxit and tol as before.
   Inside the routine an infinite right endpoint is replaced by a stand-in value beyond every finite
   endpoint, three times the largest finite endpoint plus one; a clique that reaches the stand-in is
   unbounded and is returned with right end infinity. Outputs K cliques in ql, qr, p (capacity n). */
void C_npi_turnbull(const int *n, const double *L, const double *Rin, const int *maxit,
                    const double *tol, int *K, double *ql, double *qr, double *p) {
    int nn = *n;
    double big = 0.0;
    for (int i = 0; i < nn; i++) {
        if (L[i] > big) big = L[i];
        if (isfinite(Rin[i]) && Rin[i] > big) big = Rin[i];
    }
    big = 3.0 * big + 1.0;
    double *R = (double *)malloc(sizeof(double) * nn);
    for (int i = 0; i < nn; i++) R[i] = isfinite(Rin[i]) ? Rin[i] : big;
    double *Ls = (double *)malloc(sizeof(double) * nn), *Rs = (double *)malloc(sizeof(double) * nn);
    memcpy(Ls, L, sizeof(double) * nn); memcpy(Rs, R, sizeof(double) * nn);
    int nL = uniq(Ls, nn), nR = uniq(Rs, nn), k = 0;
    for (int a = 0; a < nL; a++) for (int b = 0; b < nR; b++) {
        double l = Ls[a], r = Rs[b];
        if (r < l) continue;
        int inner = 0;
        for (int c = 0; c < nL && !inner; c++) if (Ls[c] > l && Ls[c] <= r) inner = 1;
        for (int c = 0; c < nR && !inner; c++) if (Rs[c] >= l && Rs[c] < r) inner = 1;
        if (inner) continue;
        int used = 0;                          /* keep only cliques some observation contains */
        for (int i = 0; i < nn && !used; i++) if (l >= L[i] && r <= R[i]) used = 1;
        if (used && k < nn) { ql[k] = l; qr[k] = r; k++; }
    }
    *K = k;
    if (k > 0) {
        double *pn = (double *)malloc(sizeof(double) * k), *dd = (double *)malloc(sizeof(double) * nn);
        for (int j = 0; j < k; j++) p[j] = 1.0 / k;
        int phase = 0, cap = *maxit < 2000 ? *maxit : 2000;
        for (;;) {
            for (int it = 0; it < cap; it++) {
                for (int i = 0; i < nn; i++) {
                    double s = 0.0;
                    for (int j = 0; j < k; j++) if (ql[j] >= L[i] && qr[j] <= R[i]) s += p[j];
                    dd[i] = s > 0.0 ? s : 1e-300;
                }
                double tot = 0.0;
                for (int j = 0; j < k; j++) {
                    double s = 0.0;
                    for (int i = 0; i < nn; i++) if (ql[j] >= L[i] && qr[j] <= R[i]) s += 1.0 / dd[i];
                    pn[j] = p[j] * s / nn; tot += pn[j];
                }
                double dmax = 0.0;
                for (int j = 0; j < k; j++) { pn[j] /= tot; double e = fabs(pn[j] - p[j]); if (e > dmax) dmax = e; p[j] = pn[j]; }
                if (dmax < (phase ? *tol : 1e-10)) break;
            }
            if (phase == 0 && tb_newton(nn, k, L, R, ql, qr, p)) {
                /* one self-consistency step at the optimum: no loss of precision, and the masses then satisfy the
                   update's identities exactly (for instance mass above t >= #{L_i > t} / n), as every iterate does */
                for (int i = 0; i < nn; i++) { double s = 0.0; for (int j = 0; j < k; j++) if (ql[j] >= L[i] && qr[j] <= R[i]) s += p[j]; dd[i] = s > 0.0 ? s : 1e-300; }
                double tot = 0.0;
                for (int j = 0; j < k; j++) { double s = 0.0; for (int i = 0; i < nn; i++) if (ql[j] >= L[i] && qr[j] <= R[i]) s += 1.0 / dd[i]; pn[j] = p[j] * s / nn; tot += pn[j]; }
                for (int j = 0; j < k; j++) p[j] = pn[j] / tot;
                break;
            }
            if (phase == 1) break;
            phase = 1; cap = *maxit;                                  /* not identified: the plain iteration, as before */
        }
        free(pn); free(dd);
    }
    for (int j = 0; j < k; j++) if (qr[j] >= big) qr[j] = INFINITY;   /* the stand-in was infinity */
    free(R); free(Ls); free(Rs);
}

/* ------------------------------------------------------------------------------------------- */
/* the construction                                                                             */
/* ------------------------------------------------------------------------------------------- */

/* Survival bounds of the construction at each t, with those of its members and of the Turnbull
   ambiguity set. The construction is the union of M_S and M_P; both lie inside M_C, so no clip
   is applied. Outputs, each of length nt: lo, up (construction), slo, sup (M_S), plo, pup (M_P),
   clo, cup (M_C), tlo, tup (Turnbull). Evaluate away from the data values. */
void C_npi_bounds(const int *n, const double *L, const double *R, const int *ex, const double *eps,
                  const int *closedR, const int *nt, const double *t, double *lo, double *up,
                  double *slo, double *sup, double *plo, double *pup, double *clo, double *cup,
                  double *tlo, double *tup) {
    int nn = *n, one = 1, zero = 0, maxit = 20000, K = 0;
    double tol = 1e-9, tolem = 1e-13;
    double *L2 = (double *)malloc(sizeof(double) * nn), *R2 = (double *)malloc(sizeof(double) * nn);
    double *own = (double *)malloc(sizeof(double) * (nn + 1));
    double *mlo = (double *)malloc(sizeof(double) * (2 * nn + 1)), *mhi = (double *)malloc(sizeof(double) * (2 * nn + 1));
    double *mm = (double *)malloc(sizeof(double) * (2 * nn + 1));
    double *ql = (double *)malloc(sizeof(double) * nn), *qr = (double *)malloc(sizeof(double) * nn);
    double *pp = (double *)malloc(sizeof(double) * nn);
    C_npi_untie(n, L, R, ex, eps, closedR, &tol, L2, R2);
    C_npi_share(n, L2, R2, ex, own);
    int KS = 0;
    C_npi_place(n, L2, R2, ex, own, &KS, mlo, mhi, mm);
    C_npi_mbounds(&KS, mlo, mhi, mm, &zero, nt, t, slo, sup);
    C_npi_turnbull(n, L2, R2, &maxit, &tolem, &K, ql, qr, pp);          /* the same regions for every member */
    for (int j = 0; j < *nt; j++) {
        double a = 0.0, b = 0.0;
        for (int k = 0; k < K; k++) { if (ql[k] > t[j]) a += pp[k]; if (qr[k] > t[j]) b += pp[k]; }
        tlo[j] = a; tup[j] = b;
        plo[j] = nn * a / (nn + 1.0); pup[j] = (nn * b + 1.0) / (nn + 1.0);
        int cl = 0, cu = 0;
        for (int i = 0; i < nn; i++) { if (L2[i] > t[j]) cl++; if (R2[i] > t[j]) cu++; }
        clo[j] = cl / (nn + 1.0); cup[j] = (cu + 1.0) / (nn + 1.0);
        lo[j] = slo[j] < plo[j] ? slo[j] : plo[j];
        up[j] = sup[j] > pup[j] ? sup[j] : pup[j];
    }
    (void)one;
    free(L2); free(R2); free(own); free(mlo); free(mhi); free(mm); free(ql); free(qr); free(pp);
}

/* The rc-A(n) M-functions of Coolen and Yan (2004) in closed form, for right-censored data with no
   ties; status 1 is an event, 0 a right censoring. Output n+1 intervals in time order. */
void C_npi_rc(const int *n, const double *time, const int *status, double *lo, double *hi, double *m) {
    int nn = *n, *o = (int *)malloc(sizeof(int) * (nn > 0 ? nn : 1));
    for (int i = 0; i < nn; i++) o[i] = i;
    for (int a = 1; a < nn; a++) { int v = o[a], b = a - 1; while (b >= 0 && time[o[b]] > time[v]) { o[b + 1] = o[b]; b--; } o[b + 1] = v; }
    double first = INFINITY;
    for (int i = 0; i < nn; i++) if (status[i] && time[i] < first) first = time[i];
    lo[0] = 0.0; hi[0] = first; m[0] = 1.0 / (nn + 1);
    for (int a = 0; a < nn; a++) {
        int i = o[a];
        double fac = 1.0, nx = INFINITY;
        for (int b = 0; b < a; b++) {
            int c = o[b];
            if (!status[c]) { double r = (double)(nn - b); fac *= (r + 1.0) / r; }
        }
        for (int b = a + 1; b < nn; b++) if (status[o[b]]) { nx = time[o[b]]; break; }
        lo[a + 1] = time[i]; hi[a + 1] = nx;
        m[a + 1] = status[i] ? fac / (nn + 1) : fac / ((nn + 1) * (double)(nn - a));
    }
    free(o);
}

/* ------------------------------------------------------------------------------------------- */
/* two-level NPI: exchangeable years, exchangeable units within a year                          */
/* ------------------------------------------------------------------------------------------- */

/* P(Bin(m, t) >= r), summed on the side with the smaller probability so that neither a small
   tail nor its complement is formed by cancellation. */
static double binom_tail(int m, int r, double t) {
    if (r <= 0) return 1.0;
    if (r > m) return 0.0;
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;
    double lt = log(t), l1 = log1p(-t), s = 0.0;
    if ((double)r > m * t) {
        for (int k = r; k <= m; k++) s += exp(lchoose_((double)m, (double)k) + k * lt + (m - k) * l1);
        return s > 1.0 ? 1.0 : s;
    }
    for (int k = 0; k < r; k++) s += exp(lchoose_((double)m, (double)k) + k * lt + (m - k) * l1);
    s = 1.0 - s;
    return s < 0.0 ? 0.0 : s;
}

/* Level one: past years' rates are known to lie in [L_j, R_j] within [0, 1] (L_j = R_j when the
   rate is taken as known); the construction (A(n) on exact years, the union of M_S and M_P
   otherwise) gives lower and upper survival functions S_lo, S_up for next year's rate.
   Level two: given that rate t, the number D of events among m units is Bin(m, t).
   For the increasing function g(t) = P(Bin(m, t) >= r), the lower and upper probabilities are the
   integrals of S_lo and S_up against dg over [0, 1], evaluated exactly on the segments between the
   breakpoints of the step functions (the epsilon limit of the tie convention). r vectorised. */
void C_npi_twolevel(const int *n, const double *L, const double *R, const int *m, const int *r,
                    const int *nr, const int *closedR, double *lo, double *up) {
    int nn = *n, nb = 0;
    double *bp = (double *)malloc(sizeof(double) * (2 * nn + 2));
    bp[nb++] = 0.0; bp[nb++] = 1.0;
    for (int i = 0; i < nn; i++) {
        if (L[i] > 0.0 && L[i] < 1.0) bp[nb++] = L[i];
        if (isfinite(R[i]) && R[i] > 0.0 && R[i] < 1.0) bp[nb++] = R[i];
    }
    nb = uniq(bp, nb);
    int ns = nb - 1;
    double *mid = (double *)malloc(sizeof(double) * (ns > 0 ? ns : 1));
    for (int s = 0; s < ns; s++) mid[s] = 0.5 * (bp[s] + bp[s + 1]);
    double *b[10];
    for (int k = 0; k < 10; k++) b[k] = (double *)malloc(sizeof(double) * (ns > 0 ? ns : 1));
    int *ex = (int *)malloc(sizeof(int) * (nn > 0 ? nn : 1));
    double *R1 = (double *)malloc(sizeof(double) * (nn > 0 ? nn : 1));
    for (int i = 0; i < nn; i++) { ex[i] = (L[i] == R[i]); R1[i] = isfinite(R[i]) && R[i] < 1.0 ? R[i] : 1.0; if (ex[i]) R1[i] = L[i]; }
    double eps = 0.0;                                     /* increment chosen from the data */
    C_npi_bounds(n, L, R1, ex, &eps, closedR, &ns, mid, b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9]);
    for (int j = 0; j < *nr; j++) {
        double gl = binom_tail(*m, r[j], bp[0]), a = 0.0, c = 0.0;
        for (int s = 0; s < ns; s++) {
            double gr = binom_tail(*m, r[j], bp[s + 1]), dg = gr - gl;
            a += b[0][s] * dg; c += b[1][s] * dg; gl = gr;
        }
        double g0 = binom_tail(*m, r[j], 0.0);
        lo[j] = g0 + a; up[j] = g0 + c;
        if (lo[j] > 1.0) lo[j] = 1.0;
        if (up[j] > 1.0) up[j] = 1.0;
    }
    for (int k = 0; k < 10; k++) free(b[k]);
    free(bp); free(mid); free(ex); free(R1);
}
