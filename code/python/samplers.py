"""Probability rules and the three draws (independent, EXCL, Madow).

Notation: units t = 0..n-1 in time order; `cal` is the development-data calibration (bin edges, per-bin prevalence q and variance v);
`jb` the bin index of every unit; `gaps[t]` is True when unit t starts a new recording.  Draw functions return arrays with one row per
replication: S (uploaded?), P (probability of the unit at the time of the decision) and, for the blocked designs, NEW (block starts)
and PID (period index used by the shore predictor).
"""
import math
import numpy as np

NBIN = 20
ONE16 = 65536
M32 = np.uint64(0xFFFFFFFF)
K0_PRIOR = 100.0


def calibrate(s, y, nbin=NBIN):
    """Bin the development scores into NBIN quantile bins; smoothed prevalence q and variance v per bin."""
    edges = np.quantile(s, np.linspace(0, 1, nbin + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    edges = np.unique(edges)
    k = len(edges) - 1
    idx = np.clip(np.searchsorted(edges, s, 'right') - 1, 0, k - 1)
    # smoothed bins q=(k+0.5)/(m+1) (k = sum of y in bin, m = units in bin);
    # v = mean (y-q)^2 around the smoothed q, floored at q(1-q)/(m+1) so no bin has v=0.
    q = np.zeros(k); v = np.zeros(k)
    for j in range(k):
        yj = y[idx == j]; m = len(yj)
        q[j] = (yj.sum() + 0.5) / (m + 1)
        v[j] = max(((yj - q[j]) ** 2).mean() if m else q[j] * (1 - q[j]), q[j] * (1 - q[j]) / (m + 1))
    return dict(edges=edges, q=q, v=v, qmax=float(q.max()))


def bin_of(cal, s):
    """Bin index of every score."""
    return np.clip(np.searchsorted(cal['edges'], s, 'right') - 1, 0, len(cal['q']) - 1)


def rule(q, v, mu, kappa, qmax, pmin):
    """Closed-loop rule: p = clip(sqrt(v / (mu (1 - kappa q / qmax))), pmin, 1)."""
    denom = mu * (1.0 - kappa * q / max(qmax, 1e-12))
    p = np.sqrt(v / np.maximum(denom, 1e-12))
    return np.clip(p, pmin, 1.0)


def solve_mu0(q, v, b, kappa, qmax, pmin):
    """Price mu such that the mean probability equals b on a reference score distribution."""
    lo, hi = 1e-12, 1e12
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if rule(q, v, mid, kappa, qmax, pmin).mean() > b:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def xorshift32(x):
    """Marsaglia xorshift32 (13, 17, 5), vectorised over replications; the node uses the same generator."""
    x = x ^ ((x << np.uint64(13)) & M32)
    x = x ^ (x >> np.uint64(17))
    x = x ^ ((x << np.uint64(5)) & M32)
    return x


def simulate_independent(q, v, b, kappa, qmax, pmin, mu0, eta, depth, R, rng, fixed=False, mu_levels=None):
    """Closed-loop rule with independent draws: price update and token bucket at every unit. Returns S, P, overdraft, V_plug."""
    n = len(q)
    vplug = np.zeros(R)
    S = np.zeros((R, n), bool)
    P = np.zeros((R, n))
    logmu = np.full(R, np.log(mu0))
    tokens = np.full(R, depth / 2.0)
    overdraft = np.zeros(R, bool)
    if fixed:
        state = rng.integers(1, 2 ** 32, R, dtype=np.uint64)
        lv = mu_levels  # log-mu grid for LUT
    for t in range(n):
        mu = np.exp(logmu)
        if fixed:
            k = np.clip(np.round((logmu - lv[0]) / (lv[1] - lv[0])).astype(int), 0, len(lv) - 1)
            p = rule(q[t], v[t], np.exp(lv[k]), kappa, qmax, pmin)
        else:
            p = rule(q[t], v[t], mu, kappa, qmax, pmin)
        p = np.where(tokens >= 1.0, p, pmin)          # below one token only p_min sampling
        if fixed:
            p16 = np.clip(np.round(p * 65536).astype(np.int64), 1, 65536)
            state = xorshift32(state)
            xi = (state >> np.uint64(16)).astype(np.int64) < p16   # top 16 bits vs quantised threshold
            p = p16 / 65536.0                                     # threshold actually used (uploaded)
        else:
            xi = rng.random(R) < p
        S[:, t] = xi
        P[:, t] = p
        vplug += v[t] * (1 - p) / p
        tokens = np.minimum(tokens + b - xi, depth)
        overdraft |= tokens < -1.0                   # one-sample reserve for p_min sampling
        logmu = logmu + eta * (xi - b)
    return S, P, overdraft, vplug / n ** 2


def block_flags(q, q_hot, G=2, gaps=None):
    """forced[t]: unit t must open a new block irrespective of mass (recording start, gap, or onset alignment: the bin prevalence q reaches q_hot after at least G consecutive units below it)."""
    n = len(q); hot = np.asarray(q) >= q_hot
    gaps = np.zeros(n, bool) if gaps is None else np.asarray(gaps, bool)
    forced = np.zeros(n, bool); nh = 0
    for t in range(n):
        if gaps[t]: nh = 0
        forced[t] = (t == 0) or gaps[t] or (hot[t] and nh >= G)
        nh = 0 if hot[t] else nh + 1
    return forced


def simulate_excl(q, v, q_hot, b, kappa, qmax, pmin, mu0, eta, depth, R, rng, G=2, gaps=None, fixed=False, mu_levels=None,
                  period_len=200):
    """EXCL with the closed-loop rule: blocks of mass at most 1, one uniform per block, price and token bucket updated only when a block closes. fixed=True is the integer (16-bit) path of the node."""
    n = len(q); forced = block_flags(q, q_hot, G, gaps); qm = max(qmax, 1e-12); lv = mu_levels
    S = np.zeros((R, n), bool); P = np.zeros((R, n)); NEW = np.zeros((R, n), bool); PID = np.zeros((R, n), np.int16)
    logmu = np.full(R, np.log(mu0)); tokens = np.full(R, depth / 2.0); overdraft = np.zeros(R, bool)
    plow = np.zeros(R, bool); up = np.zeros(R); blen = np.zeros(R); pid = np.zeros(R, np.int16); pstart = np.zeros(R, np.int64)
    if fixed:
        state = rng.integers(1, 2 ** 32, R, dtype=np.uint64)
        A = np.zeros(R, np.int64); u = np.zeros(R, np.int64); p16min = int(round(pmin * ONE16))
    else:
        A = np.zeros(R); u = np.ones(R)

    def mu_of(lm):
        if not fixed: return np.exp(lm)
        k = np.clip(np.round((lm - lv[0]) / (lv[1] - lv[0])).astype(int), 0, len(lv) - 1)
        return np.exp(lv[k])

    def prob(mu, plow_):
        c = 1.0 - kappa * q[t] / qm
        p = np.clip(np.sqrt(v[t] / np.maximum(mu * c, 1e-12)), pmin, 1.0)       # same as rule()
        if fixed:
            p16 = np.clip(np.round(p * ONE16).astype(np.int64), 1, ONE16)
            return np.where(plow_, p16min, p16)
        return np.where(plow_, pmin, p)

    mu_blk = mu_of(logmu)
    for t in range(n):
        p = prob(mu_blk, plow)                                          # current block's frozen state (predictable)
        new = (A + p > (ONE16 if fixed else 1.0)) | forced[t]            # boundary: depends on scores and p_t only
        if new.any():
            lm = logmu + eta * (up - b * blen); tk = np.minimum(tokens + b * blen - up, depth)
            logmu = np.where(new, lm, logmu); tokens = np.where(new, tk, tokens)
            overdraft |= new & (tokens < -1.0)
            mu_blk = np.where(new, mu_of(logmu), mu_blk); plow = np.where(new, tokens < 1.0, plow)
            p = prob(mu_blk, plow)
            A = np.where(new, 0, A); up = np.where(new, 0.0, up); blen = np.where(new, 0.0, blen)
            if fixed:
                state = np.where(new, xorshift32(state), state)
                u = np.where(new, (state >> np.uint64(16)).astype(np.int64) + 1, u)
            else:
                u = np.where(new, 1.0 - rng.random(R), u)               # (0,1]
            nst = new & (t - pstart >= period_len)
            pid = np.where(nst, pid + 1, pid).astype(np.int16); pstart = np.where(nst, t, pstart)
            NEW[:, t] = new
        An = A + p
        xi = (A < u) & (u <= An)
        up = up + xi; blen = blen + 1; A = An
        S[:, t] = xi; P[:, t] = p / ONE16 if fixed else p; PID[:, t] = pid
    return dict(S=S, P=P, overdraft=overdraft, NEW=NEW, PID=PID)


def simulate_excl_openloop(p, forced, R, rng, period_len):
    """EXCL draw at given (open-loop) probabilities p; forced marks the block starts that are not caused by the mass cap (recording starts, onset alignment)."""
    p = np.asarray(p, float); forced = np.asarray(forced, bool); n = len(p)
    assert forced[0], 'unit 0 must open a block'
    assert (p > 0).all() and (p <= 1.0).all()
    new = np.zeros(n, bool); Aprev = np.zeros(n); Acur = np.zeros(n); pidv = np.zeros(n, np.int16)
    A = 0.0; pid = 0; pstart = 0
    for t in range(n):
        pt = float(p[t])
        if forced[t] or (A + pt > 1.0):                                # same test as simulate_excl: A before adding p_t
            new[t] = True; A = 0.0
            if t - pstart >= period_len: pid += 1; pstart = t
        Aprev[t] = A; A = A + pt; Acur[t] = A; pidv[t] = pid           # sequential float addition (bit-identical)
    bid = np.cumsum(new) - 1; nb = int(bid[-1]) + 1
    U = 1.0 - rng.random((nb, R))                                      # (0,1]; rng.random((k,R)) == k successive rng.random(R) calls
    Ub = U[bid].T                                                      # (R,n): block's u at every unit
    S = (Aprev[None, :] < Ub) & (Ub <= Acur[None, :])
    return dict(S=S, P=np.broadcast_to(p, (R, n)).copy(), NEW=np.broadcast_to(new, (R, n)).copy(),
                PID=np.broadcast_to(pidv, (R, n)).copy())


def bookkeeping(p, forced, period_len):
    """Block starts and period index for probabilities p: the loop of simulate_excl_openloop."""
    p = np.asarray(p, float); forced = np.asarray(forced, bool); n = len(p)
    assert forced[0], 'unit 0 must open a block'
    assert (p > 0).all() and (p <= 1.0).all()
    new = np.zeros(n, bool); pidv = np.zeros(n, np.int16); A = 0.0; pid = 0; pstart = 0
    for t in range(n):
        pt = float(p[t])
        if forced[t] or (A + pt > 1.0):                                # same test as simulate_plugin: A before adding p_t
            new[t] = True; A = 0.0
            if t - pstart >= period_len: pid += 1; pstart = t
        A = A + pt; pidv[t] = pid
    return new, pidv


def _bcast(p, new, pid, R):
    n = len(p)
    return (np.broadcast_to(p, (R, n)).copy(), np.broadcast_to(new, (R, n)).copy(), np.broadcast_to(pid, (R, n)).copy())


def segments(gaps, n):
    seg = np.zeros(n, bool); seg[0] = True
    if gaps is not None: seg |= np.asarray(gaps, bool)
    return seg


def simulate_madow(p, gaps, forced, R, rng, period_len):
    """Madow (single start): one u in (0,1] per recording, unit t uploaded iff floor(C_t - u) > floor(C_{t-1} - u), C the cumulative probability within the recording."""
    p = np.asarray(p, float); n = len(p); seg = segments(gaps, n)
    Cp = np.zeros(n); Cc = np.zeros(n); A = 0.0
    for t in range(n):
        if seg[t]: A = 0.0                                             # C_{t-1} := 0 at the first unit of a segment
        Cp[t] = A; A = A + float(p[t]); Cc[t] = A
    sid = np.cumsum(seg) - 1; nseg = int(sid[-1]) + 1
    last = np.r_[np.flatnonzero(seg)[1:] - 1, n - 1]; segmass = Cc[last]
    U = 1.0 - rng.random((nseg, R)); Ub = U[sid].T                      # (R,n), u in (0,1]
    S = np.floor(Cc[None, :] - Ub) > np.floor(Cp[None, :] - Ub)
    new, pid = bookkeeping(p, forced, period_len); P, NEW, PID = _bcast(p, new, pid, R)
    return dict(S=S, P=P, NEW=NEW, PID=PID, SEG=seg, segmass=segmass)


def solve_scale(w_ref, b, pmin, iters=50):
    """Scale kappa with mean over the reference of clip(kappa w, pmin, 1) equal to b (bisection in log kappa)."""
    lo, hi = 1e-12, 1e12
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if np.clip(mid * w_ref, pmin, 1.0).mean() > b: hi = mid
        else: lo = mid
    return math.sqrt(lo * hi)


def causal_probs(w, w_ref, b, pmin, k0=K0_PRIOR):
    """p_t = clip(kappa* w_t m_ref / m_t, pmin, 1), m_t the running mean of past weights only (k0 prior units of weight m_ref)."""
    n = len(w); kap = solve_scale(w_ref, b, pmin); mref = float(w_ref.mean())
    cm = (k0 * mref + np.concatenate([[0.0], np.cumsum(w)[:-1]])) / (k0 + np.arange(n))
    return np.clip(kap * w * mref / cm, pmin, 1.0)


def draw_independent(p, R, rng):
    """Independent Bernoulli draws at fixed probabilities p (R replications)."""
    P = np.broadcast_to(p, (R, len(p)))
    return rng.random((R, len(p))) < P, P


def external_probs(kind, s, cal, b, ref_scores, rho=None, lam=None):
    """Open-loop probabilities from the scores only: kind 'RS' (Robust Sampling, weights sqrt(v)^(1-rho)) or 'AM' (Active Measurement, weights q)."""
    pmin = b / 4
    jb, jr = bin_of(cal, s), bin_of(cal, ref_scores)
    f, f_ref = cal['q'][jb], cal['q'][jr]; e, e_ref = np.sqrt(cal['v'][jb]), np.sqrt(cal['v'][jr])
    if kind == 'AM': return causal_probs(f, f_ref, b, pmin)
    if kind == 'RS': return causal_probs(e ** (1 - rho), e_ref ** (1 - rho), b, pmin)
    if kind == 'LAM': return causal_probs(e ** (1 - lam) * f ** lam, e_ref ** (1 - lam) * f_ref ** lam, b, pmin)
    raise ValueError(kind)
