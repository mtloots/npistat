/* npistat: nonparametric predictive inference (NPI) for censored and Bernoulli data.
 *
 * One pure-C back-end, libm only, shared byte for byte by the R package (.C interface, init.c)
 * and the Python package (ctypes). Every argument is a pointer, as the .C ABI requires.
 *
 * Two families of computation.
 *
 * 1. Censored data (Loots, "Nonparametric predictive inference for interval-censored data:
 *    M-functions anchored to the self-consistent estimator"). Each observation is known only to
 *    lie in a closed region [L_i, R_i]; R_i = +infinity for right censoring, L_i = R_i for an exact
 *    observation. Three M-function specifications:
 *      M_S  the sharing assumption ic-A(n), which is Coolen and Yan's rc-A(n) under right censoring;
 *      M_P  n/(n+1) of each Turnbull clique mass on its clique, 1/(n+1) on the whole support;
 *      M_C  the completion bound, the envelope of A(n) over all completions of the data.
 *    The construction's survival bounds are min(S_lo, P_lo) and min(max(S_up, P_up), C_up).
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

/* Hill's H(n) by the epsilon limit (Coolen and Yan 2004, Section 7): rank-ordered increments on
   the left endpoints, events before censorings at a common time, the containing region first
   among tied left endpoints; the right endpoint of a censored region is moved past an exact
   observation tied with it when closedR is 1. */
void C_npi_untie(const int *n, const double *L, const double *R, const int *ex, const double *eps,
                 const int *closedR, const double *tol, double *L2, double *R2) {
    int nn = *n, *o = (int *)malloc(sizeof(int) * (nn > 0 ? nn : 1));
    for (int i = 0; i < nn; i++) o[i] = i;
    sL_ = L; sR_ = R; sE_ = ex;
    qsort(o, nn, sizeof(int), cmp_untie);
    for (int k = 0; k < nn; k++) L2[o[k]] = L[o[k]] + (*eps) * (double)(k + 1);
    for (int i = 0; i < nn; i++) {
        if (ex[i]) { R2[i] = L2[i]; continue; }
        int tied = 0;
        if (isfinite(R[i]))
            for (int j = 0; j < nn && !tied; j++) if (ex[j] && fabs(L[j] - R[i]) < *tol) tied = 1;
        R2[i] = (*closedR && tied) ? R[i] + (*eps) * (double)(nn + 1) : R[i];
    }
    free(o);
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

/* Place the masses on intervals: observation i on (L_i, h_i), h_i the first event time after L_i,
   truncated at R_i; the region below the smallest observation on (0, first event). */
void C_npi_place(const int *n, const double *L, const double *R, const int *ex, const double *own,
                 double *lo, double *hi, double *m) {
    int nn = *n;
    double first = INFINITY;
    for (int i = 0; i < nn; i++) if (ex[i] && L[i] < first) first = L[i];
    lo[0] = 0.0; hi[0] = first; m[0] = own[0];
    for (int i = 0; i < nn; i++) {
        double nx = INFINITY;
        for (int j = 0; j < nn; j++) if (ex[j] && L[j] > L[i] && L[j] < nx) nx = L[j];
        lo[i + 1] = L[i];
        hi[i + 1] = ex[i] ? nx : (nx < R[i] ? nx : R[i]);
        m[i + 1] = own[i + 1];
    }
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

/* Maximal cliques (innermost intervals [l, r] with no other endpoint strictly inside) and the
   self-consistency iteration of Turnbull (1976). Infinite right endpoints are replaced by
   3 max(L) + 1, beyond every finite endpoint. Outputs K cliques in ql, qr, p (capacity n). */
void C_npi_turnbull(const int *n, const double *L, const double *Rin, const int *maxit,
                    const double *tol, int *K, double *ql, double *qr, double *p) {
    int nn = *n;
    double big = 0.0;
    for (int i = 0; i < nn; i++) if (L[i] > big) big = L[i];
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
        for (int it = 0; it < *maxit; it++) {
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
            if (dmax < *tol) break;
        }
        free(pn); free(dd);
    }
    free(R); free(Ls); free(Rs);
}

/* ------------------------------------------------------------------------------------------- */
/* the construction                                                                             */
/* ------------------------------------------------------------------------------------------- */

/* Survival bounds of the construction at each t, with those of its members and of the Turnbull
   ambiguity set. Outputs, each of length nt: lo, up (construction), slo, sup (M_S), plo, pup (M_P),
   clo, cup (M_C), tlo, tup (Turnbull). Evaluate away from the data values. */
void C_npi_bounds(const int *n, const double *L, const double *R, const int *ex, const double *eps,
                  const int *closedR, const int *nt, const double *t, double *lo, double *up,
                  double *slo, double *sup, double *plo, double *pup, double *clo, double *cup,
                  double *tlo, double *tup) {
    int nn = *n, one = 1, zero = 0, maxit = 20000, K = 0;
    double tol = 1e-9, tolem = 1e-13;
    double *L2 = (double *)malloc(sizeof(double) * nn), *R2 = (double *)malloc(sizeof(double) * nn);
    double *own = (double *)malloc(sizeof(double) * (nn + 1));
    double *mlo = (double *)malloc(sizeof(double) * (nn + 1)), *mhi = (double *)malloc(sizeof(double) * (nn + 1));
    double *mm = (double *)malloc(sizeof(double) * (nn + 1));
    double *ql = (double *)malloc(sizeof(double) * nn), *qr = (double *)malloc(sizeof(double) * nn);
    double *pp = (double *)malloc(sizeof(double) * nn);
    C_npi_untie(n, L, R, ex, eps, closedR, &tol, L2, R2);
    C_npi_share(n, L2, R2, ex, own);
    C_npi_place(n, L2, R2, ex, own, mlo, mhi, mm);
    int np1 = nn + 1;
    C_npi_mbounds(&np1, mlo, mhi, mm, &zero, nt, t, slo, sup);
    C_npi_turnbull(n, L, R, &maxit, &tolem, &K, ql, qr, pp);
    for (int j = 0; j < *nt; j++) {
        double a = 0.0, b = 0.0;
        for (int k = 0; k < K; k++) { if (ql[k] > t[j]) a += pp[k]; if (qr[k] > t[j]) b += pp[k]; }
        tlo[j] = a; tup[j] = b;
        plo[j] = nn * a / (nn + 1.0); pup[j] = (nn * b + 1.0) / (nn + 1.0);
        int cl = 0, cu = 0;
        for (int i = 0; i < nn; i++) { if (L[i] > t[j]) cl++; if (R[i] >= t[j]) cu++; }
        clo[j] = cl / (nn + 1.0); cup[j] = (cu + 1.0) / (nn + 1.0);
        lo[j] = slo[j] < plo[j] ? slo[j] : plo[j];
        double u = sup[j] > pup[j] ? sup[j] : pup[j];
        up[j] = u < cup[j] ? u : cup[j];
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
   rate is taken as known); the construction (A(n) on exact years, the union of M_S and M_P clipped
   to M_C otherwise) gives lower and upper survival functions S_lo, S_up for next year's rate.
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
    double eps = 1e-12;
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
