"""Exact variance tools for single-start Madow systematic sampling (open-loop plug-in).

Conventions
-----------
* Units t=0..n-1 (time ordered) are split into recording segments (contiguous units).
  ``seg_starts`` is the sorted array of first-unit indices of segments (0 always included).
* Every unit has inclusion probability p_t in (0,1].  Inside a segment the mass axis is
  C_t = sum_{s<=t, s in segment} p_s (C=0 before the segment), arc I_t = (C_{t-1}, C_t].
* Madow: one u~U(0,1] per segment; unit t is sampled iff (u+Z) meets I_t.
  On the circle R/Z the arc is (a_t, b_t] mod 1 with a_t=frac(C_{t-1}), b_t=a_t+p_t (b_t>1 => wraps).
* EXCL: segments are cut into consecutive blocks (mu_k = sum of p in block <= 1); one u_k per block,
  unit t sampled iff u_k in (A_{t-1}, A_t] (A = within-block cumulative mass).  ``blocks`` is
  the sorted array of first-unit indices of blocks (must contain every segment start).
* Estimator with fixed residuals r_t:  N(theta_hat-theta) = sum_t xi_t r_t/p_t - sum_t r_t.
  All "N^2 Var" functions below return N^2 * Var(theta_hat) (i.e. Var of sum_t xi_t r_t/p_t).

Precision: cumulative masses are accumulated in np.longdouble (80-bit, eps=1.08e-19).
``arcs`` returns a rigorous-ish worst-case bound  n_seg*eps*M_seg/2  on the error of
each segment's cumulative masses (e.g. n=2e4, M=750 -> ~8e-13 worst case; typically far smaller).
"""
import numpy as np

LD = np.longdouble
EPS_LD = float(np.finfo(LD).eps)


# ----------------------------------------------------------------------------------------
# 1. segments
# ----------------------------------------------------------------------------------------
def segments_from_gaps(n, gaps):
    """gaps: bool array length n, True = unit starts a new recording segment.
    Returns the array of segment-start indices (unit 0 is always a start)."""
    g = np.asarray(gaps, dtype=bool).copy()
    if g.shape != (n,):
        raise ValueError("gaps must have shape (n,)")
    if n > 0:
        g[0] = True
    return np.flatnonzero(g)


def _bounds(n, starts):
    starts = np.asarray(starts, dtype=np.int64)
    if len(starts) == 0 or starts[0] != 0:
        raise ValueError("starts must begin with 0")
    return starts, np.append(starts[1:], n)


def block_ids(n, starts):
    """Per-unit id of the segment/block given its start indices."""
    mark = np.zeros(n, dtype=np.int64)
    mark[np.asarray(starts, dtype=np.int64)] = 1
    return np.cumsum(mark) - 1


# ----------------------------------------------------------------------------------------
# 2. arcs
# ----------------------------------------------------------------------------------------
def arcs(p, seg_starts):
    """Arc representation of every unit on the circle R/Z.

    Returns (sid, a, b, M, errbound):
      sid[t]  segment id of unit t
      a[t]    frac(C_{t-1}) in [0,1)      (np.longdouble)
      b[t]    a[t]+p[t]  (can exceed 1 -> arc wraps around)  (np.longdouble)
      M[s]    mass of segment s            (np.longdouble)
      errbound[s]  worst-case abs error bound on cumulative masses of segment s: n_s*eps*M_s/2
    """
    p = np.asarray(p, dtype=LD)
    n = len(p)
    if np.any(p <= 0) or np.any(p > 1):
        raise ValueError("need 0<p<=1")
    st, en = _bounds(n, seg_starts)
    sid = block_ids(n, st)
    a = np.empty(n, dtype=LD)
    M = np.empty(len(st), dtype=LD)
    for s, (i, j) in enumerate(zip(st, en)):
        C = np.cumsum(p[i:j], dtype=LD)
        Cprev = np.concatenate([np.zeros(1, dtype=LD), C[:-1]])
        a[i:j] = Cprev - np.floor(Cprev)
        M[s] = C[-1]
    b = a + p
    err = (en - st).astype(float) * EPS_LD * M.astype(float) / 2
    return sid, a, b, M, err


# ----------------------------------------------------------------------------------------
# 3. exact Madow variance by boundary-event sweep
# ----------------------------------------------------------------------------------------
def var_madow(p, r, seg_starts, return_info=False):
    """Exact N^2 Var_Madow = sum_sigma Var(T_sigma),  Var(T)=int_0^1 T(u)^2 du-(int T)^2,
    T_sigma(u)=sum_{t in sigma} 1[u in (a_t,b_t] mod 1] r_t/p_t (piecewise constant in u).
    O(n log n): events (start,+w)/(end,-w) of all arcs (wrapping arcs split), sorted, cumulated.
    Returns (total, per_segment[, info]); info has max |int T - sum r| (check) and arc error bound.
    The variance is computed about the swept mean (second-moment about mean, numerically stable)."""
    p = np.asarray(p, dtype=LD)
    r = np.asarray(r, dtype=LD)
    n = len(p)
    sid, a, b, M, err = arcs(p, seg_starts)
    st, en = _bounds(n, seg_starts)
    w = r / p
    per = np.empty(len(st), dtype=LD)
    chk = 0.0
    z = np.zeros(1, dtype=LD)
    one = np.ones(1, dtype=LD)
    for s, (i, j) in enumerate(zip(st, en)):
        a_, b_, w_ = a[i:j], b[i:j], w[i:j]
        wrap = b_ > 1
        e1 = np.minimum(b_, 1)
        pos = np.concatenate([z, one, a_, e1, np.zeros(int(wrap.sum()), dtype=LD), b_[wrap] - 1])
        dl = np.concatenate([z, z, w_, -w_, w_[wrap], -w_[wrap]])
        o = np.argsort(pos, kind="stable")
        pos, dl = pos[o], dl[o]
        T = np.cumsum(dl, dtype=LD)[:-1]
        ln = np.diff(pos)
        I1 = np.sum(ln * T)
        per[s] = np.sum(ln * (T - I1) ** 2)
        chk = max(chk, abs(float(I1 - np.sum(r[i:j]))) / max(1.0, float(np.sum(np.abs(r[i:j])))))
    tot = float(np.sum(per))
    if return_info:
        return tot, per.astype(float), {"int_T_minus_sum_r_rel": chk, "arc_err_bound": err}
    return tot, per.astype(float)


# ----------------------------------------------------------------------------------------
# 4. EXCL variance, dense pairwise forms
# ----------------------------------------------------------------------------------------
def greedy_blocks(p, seg_starts, tol=1e-12):
    """Greedy deterministic blocks: start a new block at every segment start and whenever adding
    the next unit would push the block mass above 1 (+tol). Returns block-start indices."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    segset = set(int(x) for x in seg_starts)
    starts, mu = [], 0.0
    for t in range(n):
        if t in segset or mu + p[t] > 1 + tol:
            starts.append(t)
            mu = p[t]
        else:
            mu += p[t]
    return np.array(starts, dtype=np.int64)


def const_blocks(n, k, seg_starts):
    """Blocks of k consecutive units inside every segment (last block of a segment may be shorter)."""
    st, en = _bounds(n, seg_starts)
    out = []
    for i, j in zip(st, en):
        out.extend(range(i, j, k))
    return np.array(out, dtype=np.int64)


def block_masses(p, blocks):
    return np.add.reduceat(np.asarray(p, dtype=LD), np.asarray(blocks, dtype=np.int64))


def var_excl(p, r, blocks):
    """N^2 Var_EXCL = sum_t r_t^2/p_t - sum_k R_k^2 (R_k = sum of r in block k).
    (Valid whenever every block has mass <= 1.)"""
    p = np.asarray(p, dtype=LD)
    r = np.asarray(r, dtype=LD)
    R = np.add.reduceat(r, np.asarray(blocks, dtype=np.int64))
    return float(np.sum(r * r / p) - np.sum(R * R))


class _Ctx:
    pass


def _madow_ctx(p, seg_starts, dtype=LD):
    c = _Ctx()
    c.n = len(p)
    c.p = np.asarray(p, dtype=LD).astype(dtype)
    sid, a, b, M, err = arcs(p, seg_starts)
    c.sid = sid
    a = a.astype(dtype)
    b = b.astype(dtype)
    zero = np.zeros_like(a)
    c.pieces = [(a, np.minimum(b, 1)), (zero, np.where(b > 1, b - 1, zero))]
    return c


def rho_madow_rows(ctx, i0, i1):
    """Rows i0..i1 of rho_st = pi_st/(p_s p_t) for Madow (same segment: arc overlap mod 1;
    different segments: 1; diagonal 1/p_t)."""
    ov = 0
    for (sa, ea) in ctx.pieces:
        for (sb, eb) in ctx.pieces:
            ov = ov + np.maximum(0, np.minimum(ea[i0:i1, None], eb[None, :])
                                 - np.maximum(sa[i0:i1, None], sb[None, :]))
    rho = ov / (ctx.p[i0:i1, None] * ctx.p[None, :])
    rho = np.where(ctx.sid[i0:i1, None] == ctx.sid[None, :], rho, 1)
    idx = np.arange(i0, i1)
    rho[idx - i0, idx] = 1 / ctx.p[idx]
    return rho


def _excl_ctx(p, blocks, dtype=LD):
    c = _Ctx()
    c.n = len(p)
    c.p = np.asarray(p, dtype=LD).astype(dtype)
    c.bid = block_ids(len(p), blocks)
    return c


def rho_excl_rows(ctx, i0, i1):
    """Rows of rho for EXCL: same block 0, different block 1, diagonal 1/p."""
    same = ctx.bid[i0:i1, None] == ctx.bid[None, :]
    rho = np.where(same, 0, 1).astype(ctx.p.dtype)
    idx = np.arange(i0, i1)
    rho[idx - i0, idx] = 1 / ctx.p[idx]
    return rho


def rho_madow(p, seg_starts, dtype=LD):
    """Dense rho matrix of Madow (n<=3000 advisable)."""
    c = _madow_ctx(p, seg_starts, dtype)
    return rho_madow_rows(c, 0, c.n)


def rho_excl(p, blocks, dtype=LD):
    c = _excl_ctx(p, blocks, dtype)
    return rho_excl_rows(c, 0, c.n)


def pairwise_var(rho, r):
    """N^2 Var = sum_{s,t} r_s r_t (rho_st - 1); rho_st=pi_st/(p_s p_t) (s!=t), rho_tt=1/p_t."""
    rho = np.asarray(rho)
    r = np.asarray(r, dtype=rho.dtype)
    return float(r @ ((rho - 1) @ r))


def pairwise_chunked(rho_rows, n, weight_rows, chunk=200):
    """sum_{s,t} W_st (rho_st - 1) in row chunks; rho_rows(i0,i1), weight_rows(i0,i1) -> (i1-i0,n)."""
    tot = 0
    for i0 in range(0, n, chunk):
        i1 = min(n, i0 + chunk)
        tot = tot + np.sum(weight_rows(i0, i1) * (rho_rows(i0, i1) - 1))
    return float(tot)


def expected_var(p, seg_starts, acov, kind="madow", blocks=None, chunk=200, dtype=LD):
    """E_zeta N^2 Var for random residuals with E r_s r_t = acov[|s-t|]:
    sum_{s,t} (rho_st-1) acov[|s-t|].  kind in {'madow','excl'}."""
    n = len(p)
    acov = np.asarray(acov, dtype=dtype)
    if kind == "madow":
        ctx = _madow_ctx(p, seg_starts, dtype)
        rows = lambda i0, i1: rho_madow_rows(ctx, i0, i1)
    else:
        ctx = _excl_ctx(p, blocks, dtype)
        rows = lambda i0, i1: rho_excl_rows(ctx, i0, i1)
    jj = np.arange(n)
    wr = lambda i0, i1: acov[np.abs(np.arange(i0, i1)[:, None] - jj[None, :])]
    return pairwise_chunked(rows, n, wr, chunk)


# ----------------------------------------------------------------------------------------
# 5. covariogram terms
# ----------------------------------------------------------------------------------------
def _gamma_h(C, g, h):
    """gamma(h)=int g(x)g(x+h)dx for piecewise-constant g (value g[k] on (C[k],C[k+1]]); exact."""
    M = C[-1]
    if h >= M:
        return LD(0)
    top = M - h
    pts = np.concatenate([np.zeros(1, dtype=LD), np.array([top], dtype=LD),
                          C[C <= top], (C - h)[(C - h) >= 0]])
    pts = np.unique(pts)
    pts = pts[(pts >= 0) & (pts <= top)]
    ln = np.diff(pts)
    mid = (pts[:-1] + pts[1:]) / 2
    m = len(g)
    i1 = np.clip(np.searchsorted(C, mid, side="left") - 1, 0, m - 1)
    i2 = np.clip(np.searchsorted(C, mid + h, side="left") - 1, 0, m - 1)
    return np.sum(ln * g[i1] * g[i2])


def covariogram_terms(p, r, seg_starts, jmax=5):
    """Covariogram decomposition  N^2 Var_M = sum_sigma gamma_sigma(0) - Delta  with
       gamma_sigma(0)=sum r^2/p,  int_0^inf gamma_sigma=(sum_sigma r)^2/2,
       Delta = sum_sigma 2[int_0^inf gamma_sigma - sum_{j>=1} gamma_sigma(j)]
             = sum_sigma[(sum r)^2 - 2 sum_{j>=1} gamma_sigma(j)].
    Returns dict: gamma0, gamma0_seg, gamma_j (sum over segments, j=1..jmax), gamma_sum_seg
    (sum_{j>=1} gamma(j), all j<M), int_pos_seg, Delta_seg, Delta, var (=gamma0-Delta), var_seg."""
    p = np.asarray(p, dtype=LD)
    r = np.asarray(r, dtype=LD)
    n = len(p)
    st, en = _bounds(n, seg_starts)
    g0 = np.empty(len(st), dtype=LD)
    gs = np.empty(len(st), dtype=LD)
    ip = np.empty(len(st), dtype=LD)
    gj = np.zeros(jmax, dtype=LD)
    for s, (i, j) in enumerate(zip(st, en)):
        pp, rr = p[i:j], r[i:j]
        g = rr / pp
        C = np.concatenate([np.zeros(1, dtype=LD), np.cumsum(pp, dtype=LD)])
        g0[s] = np.sum(rr * rr / pp)
        ip[s] = np.sum(rr) ** 2 / 2
        tot = LD(0)
        for jj in range(1, int(np.floor(C[-1])) + 1):
            v = _gamma_h(C, g, LD(jj))
            tot += v
            if jj <= jmax:
                gj[jj - 1] += v
        gs[s] = tot
    D_seg = 2 * (ip - gs)
    D = float(np.sum(D_seg))
    return {"gamma0": float(np.sum(g0)), "gamma0_seg": g0.astype(float), "gamma_j": gj.astype(float),
            "gamma_sum_seg": gs.astype(float), "int_pos_seg": ip.astype(float),
            "Delta_seg": D_seg.astype(float), "Delta": D,
            "var": float(np.sum(g0)) - D, "var_seg": (g0 - D_seg).astype(float)}


# ----------------------------------------------------------------------------------------
# 6. Gamma
# ----------------------------------------------------------------------------------------
def gamma_cross(p, r, seg_starts, blocks):
    """Gamma = N^2 [Var_M - Var_E]  (exact; uses var_madow and var_excl)."""
    return var_madow(p, r, seg_starts)[0] - var_excl(p, r, blocks)
