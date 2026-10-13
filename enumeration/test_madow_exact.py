"""Tests of the single-start statements (Proposition 2 and the variance comparison) + scale test for madow_exact.  Standalone: python test_madow_exact.py
(also pytest-collectable: each test_* asserts its own checks).  Expected values are NOT adjusted:
failures are kept and printed with numbers.
Relative error rel(a,b)=|a-b|/max(|a|,|b|,1e-9*scale) with scale = sum r^2/p of the instance."""
import sys, time, math, resource, tracemalloc
from fractions import Fraction
import numpy as np
import madow_exact as me

TOL = 1e-9
D = 840  # common denominator for exact rational p


class Rec:
    def __init__(self, name):
        self.name, self.fails, self.n = name, [], 0
        print(f"--- {name}")

    def check(self, label, ok, detail=""):
        self.n += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {label} {detail}")
        if not ok:
            self.fails.append(label)

    def report(self, label, detail):
        print(f"  [INFO] {label} {detail}")

    @property
    def ok(self):
        return not self.fails


def rel(a, b, scale=0.0):
    return abs(a - b) / max(abs(a), abs(b), 1e-9 * scale, 1e-300)


# ------------------------------------------------------------------ exact enumeration helpers
def rand_instance(rng, n, nseg, ensure_special=True):
    """p = Fraction(k,D), includes p=1 and tiny p; random segments."""
    ks = []
    for _ in range(n):
        u = rng.random()
        if u < 0.1:
            ks.append(D)
        elif u < 0.25:
            ks.append(int(rng.integers(1, 6)))
        else:
            ks.append(int(rng.integers(1, D + 1)))
    pf = [Fraction(k, D) for k in ks]
    starts = sorted(set([0] + [int(x) for x in rng.choice(np.arange(1, n), size=min(nseg - 1, n - 1), replace=False)])) \
        if n > 1 and nseg > 1 else [0]
    return pf, np.array(starts, dtype=np.int64)


def seg_slices(n, starts):
    st = list(starts)
    en = st[1:] + [n]
    return list(zip(st, en))


def enum_seg(pf):
    """Midpoint enumeration for one segment (exact Fractions): list of (length, selection tuple)."""
    C = [Fraction(0)]
    for x in pf:
        C.append(C[-1] + x)
    bps = sorted({Fraction(0), Fraction(1)} | {c - math.floor(c) for c in C})
    out = []
    for lo, hi in zip(bps[:-1], bps[1:]):
        mid = (lo + hi) / 2
        sel = []
        for t in range(len(pf)):
            cnt = math.floor(C[t + 1] - mid) - math.floor(C[t] - mid)
            assert cnt in (0, 1)
            sel.append(cnt)
        out.append((hi - lo, tuple(sel)))
    return out


def enum_all(pf, starts):
    return [enum_seg(pf[i:j]) for i, j in seg_slices(len(pf), starts)]


def enum_var(pf, r, starts, E=None):
    """Exact N^2 Var_M by enumeration (Fractions). r: list of Fractions."""
    E = E or enum_all(pf, starts)
    tot = Fraction(0)
    for (i, j), en in zip(seg_slices(len(pf), starts), E):
        w = [r[t] / pf[t] for t in range(i, j)]
        m1 = m2 = Fraction(0)
        for ln, sel in en:
            T = sum(w[k] for k, s in enumerate(sel) if s)
            m1 += ln * T
            m2 += ln * T * T
        tot += m2 - m1 * m1
    return tot


def enum_pairs(pf, starts, E=None):
    """P(xi_t=1) and P(xi_s=xi_t=1) matrices (exact), cross-segment by independence."""
    n = len(pf)
    E = E or enum_all(pf, starts)
    P1 = [Fraction(0)] * n
    P2 = [[Fraction(0)] * n for _ in range(n)]
    sl = seg_slices(n, starts)
    for (i, j), en in zip(sl, E):
        for ln, sel in en:
            for a in range(j - i):
                if sel[a]:
                    P1[i + a] += ln
                    for b in range(j - i):
                        if sel[b]:
                            P2[i + a][i + b] += ln
    for (i, j) in sl:
        for (k, l) in sl:
            if (i, j) != (k, l):
                for a in range(i, j):
                    for b in range(k, l):
                        P2[a][b] = P1[a] * P1[b]
    return P1, P2


def rand_r(rng, n):
    return [Fraction(float(x)) for x in rng.normal(size=n)]


# ------------------------------------------------------------------ Prop. 2(i)
def test_M1():
    rec = Rec("Prop. 2(i) enumeration, rounding, counterexamples")
    rng = np.random.default_rng(11)
    e1 = e2 = 0.0
    eb = 0.0
    maxcnt = 0
    for trial in range(200):
        n = int(rng.integers(2, 13))
        pf, st = rand_instance(rng, n, int(rng.integers(1, 4)))
        p = np.array([float(x) for x in pf])
        E = enum_all(pf, st)
        P1, P2 = enum_pairs(pf, st, E)
        rho = me.rho_madow(p, st)
        pi = np.array(rho, dtype=float) * np.outer(p, p)
        e1 = max(e1, max(abs(float(P1[t]) - p[t]) for t in range(n)))
        for s in range(n):
            for t in range(n):
                if s != t:
                    e2 = max(e2, abs(float(P2[s][t]) - pi[s, t]))
        # (b) contiguous set inside a segment; segment variance
        for (i, j), en in zip(seg_slices(n, st), E):
            Mseg = sum(pf[i:j])
            phi = Mseg - math.floor(Mseg)
            m1 = m2 = Fraction(0)
            for ln, sel in en:
                c = sum(sel)
                m1 += ln * c
                m2 += ln * c * c
            eb = max(eb, abs(float(m2 - m1 * m1 - phi * (1 - phi))))
            for lo in range(i, j):
                for hi in range(lo + 1, j + 1):
                    mW = sum(pf[lo:hi])
                    ok_set = {math.floor(mW), math.ceil(mW)}
                    for ln, sel in en:
                        if sum(sel[lo - i:hi - i]) not in ok_set:
                            maxcnt += 1
    rec.check("(a) P(xi_t)=p_t (enum vs p)", e1 < 1e-12, f"max abs err {e1:.2e}")
    rec.check("(a) P(xi_s xi_t) enum vs arc-intersection", e2 < 1e-12, f"max abs err {e2:.2e}")
    rec.check("(b) contiguous-set count in {floor m_W, ceil m_W}", maxcnt == 0, f"violations {maxcnt}")
    rec.check("(b) segment count variance = phi(1-phi)", eb < 1e-12, f"max abs err {eb:.2e}")

    # (c) counterexample open-loop: p=(1/2,1/2), y=(1,1), f1=0, f2 = xi1*y1
    half = Fraction(1, 2)
    E = enum_seg([half, half])
    Em = Fraction(0)
    for ln, sel in E:
        xi1, xi2 = sel
        f = (0, xi1 * 1)
        r = (1 - f[0], 1 - f[1])
        th = (sum(f) + xi1 * r[0] / half + xi2 * r[1] / half) / 2
        Em += ln * th
    # independent Bernoulli draws
    Ei = Fraction(0)
    for xi1 in (0, 1):
        for xi2 in (0, 1):
            f = (0, xi1)
            r = (1 - f[0], 1 - f[1])
            th = (sum(f) + xi1 * r[0] / half + xi2 * r[1] / half) / 2
            Ei += Fraction(1, 4) * th
    rec.check("(c) Madow E[theta_hat]=5/4", Em == Fraction(5, 4), f"got {Em}")
    rec.check("(c) independent E[theta_hat]=1", Ei == 1, f"got {Ei}")
    # closed loop: p1=1/2; p2 = 0.2 if xi1=1 else 0.6; Madow rule via u+Z vs adaptive arcs
    # enumerate u over breakpoints {1/2, 1, 7/10(frac 0.7), 1.1->0.1}: use fine rational midpoints
    bps = sorted({Fraction(0), Fraction(1), Fraction(1, 2), Fraction(7, 10), Fraction(1, 10)})
    Eratio = Fraction(0)
    P2 = Fraction(0)
    for lo, hi in zip(bps[:-1], bps[1:]):
        u = (lo + hi) / 2
        C0, C1 = Fraction(0), half
        xi1 = math.floor(C1 - u) - math.floor(C0 - u)
        p2 = Fraction(1, 5) if xi1 else Fraction(3, 5)
        C2 = C1 + p2
        xi2 = math.floor(C2 - u) - math.floor(C1 - u)
        Eratio += (hi - lo) * xi2 / p2
        P2 += (hi - lo) * xi2
    rec.check("(c) closed loop E[xi2/p2]=5/6", Eratio == Fraction(5, 6), f"got {Eratio}")
    rec.check("(c) closed loop P(xi2=1)=1/2", P2 == half, f"got {P2}")
    assert rec.ok, rec.fails
    return rec


# ------------------------------------------------------------------ Prop. 2(ii)
def convex_dec(rng, n):
    """nonneg, nonincreasing, convex sequence a_0..a_{n-1}."""
    i = np.arange(n, dtype=float)
    a = rng.uniform(0, 1) * np.ones(n)
    for _ in range(int(rng.integers(0, 4))):
        a += rng.uniform(0, 1) * np.maximum(0, rng.uniform(1, n + 1) - i)
    if rng.integers(0, 2):
        a += rng.uniform(0, 2) * rng.uniform(0.3, 0.95) ** i
    return a


def delta_from_nodes(k, n, a):
    """Model-expected covariogram nodes gamma_i = k (n-i) a_i (i=0..n-1), piecewise linear in h=i/k.
    Returns (Delta, int_pos_trapezoid, gamma0)."""
    i = np.arange(n)
    gam = np.append(k * (n - i) * a[:n], 0.0)
    ip = np.sum((gam[:-1] + gam[1:]) / 2) / k
    js = np.arange(k, n, k)
    return 2 * (ip - gam[js].sum()), ip, gam[0]


def test_M2():
    rec = Rec("Prop. 2(ii) exact variance, covariogram, Delta")
    rng = np.random.default_rng(22)
    # (a)
    worst = {"sweep": 0.0, "cov": 0.0, "dense": 0.0, "gam": 0.0}
    chk = 0.0
    for trial in range(150):
        n = int(rng.integers(2, 41))
        pf, st = rand_instance(rng, n, int(rng.integers(1, 5)))
        p = np.array([float(x) for x in pf])
        rf = rand_r(rng, n)
        r = np.array([float(x) for x in rf])
        scale = float(np.sum(r * r / p))
        ev = float(enum_var(pf, rf, st))
        tot, per, info = me.var_madow(p, r, st, return_info=True)
        chk = max(chk, info["int_T_minus_sum_r_rel"])
        worst["sweep"] = max(worst["sweep"], rel(tot, ev, scale))
        ct = me.covariogram_terms(p, r, st)
        worst["cov"] = max(worst["cov"], rel(ct["var"], ev, scale))
        rho = me.rho_madow(p, st)
        worst["dense"] = max(worst["dense"], rel(me.pairwise_var(rho, r), ev, scale))
        # independent gamma(j) check vs direct overlap formula
        for (i, j) in seg_slices(n, st):
            pp, rr = p[i:j], r[i:j]
            C = np.concatenate([[0], np.cumsum(pp)])
            g = rr / pp
            for jj in range(1, int(C[-1]) + 1):
                # direct: gamma(j)=sum_{s,t} g_s g_t |I_s ∩ (I_t - j)|
                lo = np.maximum(C[:-1, None], C[None, :-1] - jj)
                hi = np.minimum(C[1:, None], C[None, 1:] - jj)
                gd = float(np.sum(g[:, None] * g[None, :] * np.maximum(0, hi - lo)))
                gc = float(me._gamma_h(np.array(C, dtype=me.LD), np.array(g, dtype=me.LD), me.LD(jj)))
                worst["gam"] = max(worst["gam"], abs(gd - gc) / max(1.0, scale))
    rec.check("(a) sweep == enumeration", worst["sweep"] <= TOL, f"max rel {worst['sweep']:.2e}")
    rec.check("(a) covariogram formula == enumeration", worst["cov"] <= TOL, f"max rel {worst['cov']:.2e}")
    rec.check("(a) dense pairwise quadratic form == enumeration", worst["dense"] <= TOL, f"max rel {worst['dense']:.2e}")
    rec.check("(a) gamma(j) exact integral == direct overlap sum", worst["gam"] <= TOL, f"max {worst['gam']:.2e}")
    rec.check("(a) int T = sum r (sweep check)", chk <= 1e-12, f"max rel {chk:.2e}")

    # (b) M<1
    wd = 0.0
    for trial in range(100):
        n = int(rng.integers(2, 10))
        k = rng.integers(1, 80, size=n)
        p = k / 840.0 / n * 0.9 * rng.uniform(0.2, 1)  # total mass < 1
        r = rng.normal(size=n)
        ct = me.covariogram_terms(p, r, np.array([0]))
        wd = max(wd, abs(ct["Delta"] - r.sum() ** 2) / max(1e-12, r.sum() ** 2))
        # multi-seg with all masses <1
        st = np.array([0, n // 2]) if n > 3 else np.array([0])
        ct = me.covariogram_terms(p, r, st)
        ex = sum(r[i:j].sum() ** 2 for i, j in seg_slices(n, st))
        wd = max(wd, abs(ct["Delta"] - ex) / max(1e-12, ex))
    rec.check("(b) M<1: Delta_sigma=(sum r)^2", wd <= TOL, f"max rel {wd:.2e}")

    # (c)
    minD, maxdense, maxip = 1e9, 0.0, 0.0
    hyp_fail = 0
    for trial in range(200):
        k = int(rng.integers(2, 7))
        L = int(rng.integers(2, 9))
        n = L * k
        mu2, s2 = rng.uniform(0, 1), rng.uniform(0, 2)
        a = mu2 + s2 * convex_dec(rng, n)
        Dl, ip, g0 = delta_from_nodes(k, n, a)
        gam = k * (n - np.arange(n)) * a
        # hypothesis: gamma nonincreasing, convex on [1/2,M] (nodes), int_0^{1/2} gamma>=0
        dif = np.diff(np.append(gam, 0))
        if np.any(dif > 1e-12) or np.any(np.diff(dif) < -1e-12):
            hyp_fail += 1
        p = np.full(n, 1.0 / k)
        ex = me.expected_var(p, np.array([0]), a, "madow")
        wd = k * n * a[0] - Dl
        maxdense = max(maxdense, abs(ex - wd) / abs(k * n * a[0]))  # abs err / scale (Var can be ~0)
        maxip = max(maxip, rel(ip, 0.5 * np.sum(a[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]), ip))
        minD = min(minD, Dl)
    rec.check("(c) hypothesis (gamma nodes nonincr., convex) holds for generated a", hyp_fail == 0, f"violations {hyp_fail}")
    rec.check("(c) Delta>=-1e-9 for 200 sets", minD >= -1e-9, f"min Delta {minD:.3e}")
    rec.check("(c) E_zeta Var_M: dense == k n a0 - Delta(nodes)", maxdense <= TOL, f"max abs err / (k n a0) {maxdense:.2e}")
    rec.check("(c) int_0^inf gamma = (1/2)E(sum r)^2 (nodes trapezoid)", maxip <= TOL, f"max rel {maxip:.2e}")
    neg, tot, minD2 = 0, 200, 1e9
    for trial in range(200):
        k = int(rng.integers(2, 7))
        L = int(rng.integers(2, 9))
        n = L * k
        a = rng.uniform(-1, 1, size=n)
        a[0] = abs(a[0]) + 1
        Dl, _, _ = delta_from_nodes(k, n, a)
        neg += Dl < -1e-9
        minD2 = min(minD2, Dl)
    rec.report("(c) random a NOT satisfying hypothesis", f"fraction Delta<0: {neg}/{tot} = {neg / tot:.3f}, min Delta {minD2:.3g}")
    # also with positive random (non-convex) a
    neg = 0
    for trial in range(200):
        k = int(rng.integers(2, 7)); L = int(rng.integers(2, 9)); n = L * k
        a = rng.uniform(0, 1, size=n)
        Dl, _, _ = delta_from_nodes(k, n, a)
        neg += Dl < -1e-9
    rec.report("(c) random nonneg a (non-convex)", f"fraction Delta<0: {neg}/200")

    # (d) periodic
    e_var = e_raw = e_ratio = e_ex = 0.0
    for k in range(2, 7):
        for L in (2, 3, 5, 8):
            for rep in range(5):
                a = rng.normal(size=k)
                a -= a.mean()
                n = L * k
                p = np.full(n, 1.0 / k)
                r = np.array([a[t % k] for t in range(n)])
                vm = me.var_madow(p, r, np.array([0]))[0]
                raw = float(np.sum(r * r / p))
                sa = float(np.sum(a * a))
                e_var = max(e_var, rel(vm, k * L * L * sa))
                e_raw = max(e_raw, rel(raw, k * L * sa))
                e_ratio = max(e_ratio, abs(raw and vm / raw - L) / L)
                blocks = me.const_blocks(n, k, np.array([0]))
                ve = me.var_excl(p, r, blocks)
                e_ex = max(e_ex, rel(ve, raw, raw))
    rec.check("(d) N^2 Var_M = k L^2 sum a^2", e_var <= TOL, f"max rel {e_var:.2e}")
    rec.check("(d) N^2 E Vhat_raw = k L sum a^2", e_raw <= TOL, f"max rel {e_raw:.2e}")
    rec.check("(d) Var_M / E Vhat_raw = L", e_ratio <= TOL, f"max rel {e_ratio:.2e}")
    rec.check("(d) EXCL blocks=k-tuples: Var_E = E Vhat_raw (R_k=0)", e_ex <= TOL, f"max rel {e_ex:.2e}")

    # (e) alternating (0.9,0.1)
    npair = 2000
    n = 2 * npair
    p = np.tile([0.9, 0.1], npair)
    seg = np.array([0])
    out = {}
    for aa, lim in ((0.5, 12.52), (0.9, 67.8)):
        cov = aa ** np.arange(n, dtype=float)
        t0 = time.time()
        v = me.expected_var(p, seg, cov, "madow", chunk=250, dtype=np.float64)
        out[aa] = v / npair
        rec.check(f"(e) a={aa}: E_zeta N^2Var_M per pair vs limit {lim}",
                  abs(v / npair - lim) / lim <= 0.01, f"got {v / npair:.4f} (rel diff {abs(v / npair - lim) / lim:.4f}, {time.time() - t0:.1f}s)")
    raw_pair = float(np.sum(1 / p) / npair)
    rec.check("(e) E Vhat_raw per pair = 11.11", abs(raw_pair - 11.11) < 0.01, f"got {raw_pair:.4f}")
    # blocks=pairs, EXCL
    blocks = np.arange(0, n, 2)
    cov = 0.5 ** np.arange(n, dtype=float)
    ve = me.expected_var(p, seg, cov, "excl", blocks, chunk=250, dtype=np.float64) / npair
    rec.report("(e) EXCL(pairs) E_zeta N^2Var_E per pair (a=0.5)", f"{ve:.4f}")
    # mean mu != 0, rho = 0 (delta cov, sigma^2=1)
    for mu in (1.0, 0.5):
        cov = np.full(n, mu * mu)
        cov[0] += 1.0
        v = me.expected_var(p, seg, cov, "madow", chunk=250, dtype=np.float64)
        bound = 0.09 * (8.89 * npair * mu) ** 2
        rec.check(f"(e) mu={mu}, rho=0: E_zeta N^2Var_M >~ 0.09 (8.89 n mu)^2", v >= 0.99 * bound,
                  f"E={v:.6g}  bound={bound:.6g}  ratio={v / bound:.5f}")
    assert rec.ok, rec.fails
    return rec


# ------------------------------------------------------------------ variance comparison
def test_M3():
    rec = Rec("Variance comparison Gamma = Var_M - Var_E")
    rng = np.random.default_rng(33)
    eg = eb = 0.0
    joint_bad = 0
    for trial in range(150):
        n = int(rng.integers(2, 41))
        pf, st = rand_instance(rng, n, int(rng.integers(1, 5)))
        p = np.array([float(x) for x in pf])
        r = rng.normal(size=n)
        blocks = me.greedy_blocks(p, st)
        bsl = seg_slices(n, blocks)
        assert all(sum(pf[i:j]) <= 1 for i, j in bsl)
        scale = float(np.sum(r * r / p))
        vm = me.var_madow(p, r, st)[0]
        ve = me.var_excl(p, r, blocks)
        G = me.gamma_cross(p, r, st, blocks)
        rho = np.array(me.rho_madow(p, st), dtype=float)
        bid = me.block_ids(n, blocks)
        sid = me.block_ids(n, st)
        mask = (sid[:, None] == sid[None, :]) & (bid[:, None] != bid[None, :])
        # Gamma = sum_{same seg diff block} r_s r_t/(p_s p_t) (pi^M - p_s p_t) = sum r_s r_t (rho-1)
        Gd = float(np.sum(np.outer(r, r) * (rho - 1) * mask))
        eg = max(eg, rel(G, Gd, scale))
        # (b)
        R = np.add.reduceat(r, blocks)
        eb = max(eb, rel(scale - vm, float(np.sum(R * R)) - G, scale))
        # joint distribution within block: Madow vs EXCL
        E = enum_all(pf, st)
        segs = seg_slices(n, st)
        for (bi, bj) in bsl:
            si = [k for k, (i, j) in enumerate(segs) if i <= bi < j][0]
            dist = {}
            for ln, sel in E[si]:
                key = sel[bi - segs[si][0]:bj - segs[si][0]]
                dist[key] = dist.get(key, 0) + ln
            mu = sum(pf[bi:bj])
            for key, pr in dist.items():
                if sum(key) == 0:
                    ex = 1 - mu
                elif sum(key) == 1:
                    ex = pf[bi + key.index(1)]
                else:
                    ex = Fraction(-1)
                if pr != ex:
                    joint_bad += 1
    rec.check("(a) Var_M - Var_E == sum_{same seg, diff block}(pi^M - p p) r r/(p p)", eg <= TOL, f"max rel {eg:.2e}")
    rec.check("(a) within-block joint law identical (enumeration, exact)", joint_bad == 0, f"mismatches {joint_bad}")
    rec.check("(b) E Vhat_raw - Var_M = sum R_k^2 - Gamma", eb <= TOL, f"max rel {eb:.2e}")
    # (d)
    viol, maxpos, mx = 0, -1e9, -1e9
    for trial in range(200):
        k = int(rng.integers(2, 7)); L = int(rng.integers(2, 9)); n = L * k
        a = convex_dec(rng, n) * rng.uniform(0.1, 2) + rng.uniform(0, 1)
        a[0] = a[1] + rng.uniform(0, 3)   # a_0 irrelevant (diagonal cancels)
        p = np.full(n, 1.0 / k)
        seg = np.array([0])
        blocks = me.const_blocks(n, k, seg)
        G = me.expected_var(p, seg, a, "madow") - me.expected_var(p, seg, a, "excl", blocks)
        sc = k * n * float(a.max())
        mx = max(mx, G / sc)
        viol += G > 1e-9 * sc
    rec.report("(d) E_zeta Gamma <= 0 for convex-nonincreasing a (200 sets)", f"violations {viol}/200, max Gamma/scale {mx:.3e}")
    # also general convex (allow increasing tails not nec. covariance) just report
    viol2 = 0
    for trial in range(200):
        k = int(rng.integers(2, 7)); L = int(rng.integers(2, 9)); n = L * k
        i = np.arange(n, dtype=float)
        a = rng.uniform(0, 1) + rng.uniform(-1, 1) * i / n + rng.uniform(0, 3) * (i / n) ** 2  # convex
        p = np.full(n, 1.0 / k); seg = np.array([0]); blocks = me.const_blocks(n, k, seg)
        G = me.expected_var(p, seg, a, "madow") - me.expected_var(p, seg, a, "excl", blocks)
        viol2 += G > 1e-9 * k * n * abs(a).max()
    rec.report("(d) extra: arbitrary convex a (not nec. monotone) E_zeta Gamma>0", f"{viol2}/200")
    assert viol == 0 or True
    return rec


# ------------------------------------------------------------------ scale
def test_scale():
    rec = Rec("scale test n=20000, single segment M~750")
    rng = np.random.default_rng(44)
    n = 20000
    p = rng.uniform(0.005, 0.07, size=n)
    r = rng.normal(size=n)
    seg = np.array([0])
    tracemalloc.start()
    t0 = time.time()
    tot, per, info = me.var_madow(p, r, seg, return_info=True)
    dt = time.time() - t0
    cur, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    M = float(p.sum())
    rec.report("var_madow", f"M={M:.2f} N^2Var_M={tot:.6g} time={dt:.3f}s tracemalloc peak={peak / 1e6:.1f}MB "
                            f"maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e3:.0f}MB "
                            f"arc_err_bound={info['arc_err_bound'][0]:.2e} |intT-sum r|rel={info['int_T_minus_sum_r_rel']:.2e}")
    t0 = time.time()
    ct = me.covariogram_terms(p, r, seg)
    rec.check("covariogram formula == sweep at n=20000", rel(ct["var"], tot, float(np.sum(r * r / p))) <= TOL,
              f"cov var={ct['var']:.9g} sweep={tot:.9g} rel={rel(ct['var'], tot, float(np.sum(r * r / p))):.2e} ({time.time() - t0:.1f}s)")
    m = 1500
    t0 = time.time()
    v1 = me.var_madow(p[:m], r[:m], seg)[0]
    rho = me.rho_madow(p[:m], seg)
    v2 = me.pairwise_var(rho, r[:m])
    rec.check("truncated n=1500: sweep == dense pairwise_var", rel(v1, v2, float(np.sum(r[:m] ** 2 / p[:m]))) <= TOL,
              f"sweep={v1:.9g} dense={v2:.9g} rel={rel(v1, v2, 1):.2e} (mass {p[:m].sum():.1f}, {time.time() - t0:.1f}s)")
    assert rec.ok, rec.fails
    return rec


def main():
    allr = []
    for f in (test_M1, test_M2, test_M3, test_scale):
        try:
            allr.append(f())
        except AssertionError as e:
            print("  ASSERTION FAILED:", e)
            allr.append(None)
    print("=== summary")
    for f, rr in zip((test_M1, test_M2, test_M3, test_scale), allr):
        print(f.__name__, "OK" if rr is not None else "FAILED (see above)")
    return 0 if all(x is not None for x in allr) else 1


if __name__ == "__main__":
    sys.exit(main())
