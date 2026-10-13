"""Shore estimator, variance estimates and interval helpers.

theta_hat = N^-1 sum_t [f_t + xi_t r_t / p_t], r_t = y_t - f_t, with the prequential predictor f_t = sigmoid(logit q_j(t) + a_k):
the intercept a_k of period k is a ridge-penalised logistic fit to the weighted uploads of earlier periods.
"""
import numpy as np

Z95 = 1.959964


def _sig(x): return 1.0 / (1.0 + np.exp(-x))


def _logit(q): return np.log(q) - np.log1p(-q)


def fit_intercept(lq, W1, W0, a0, shrink, iters=12):
    """argmax_a sum_j W1 log q_j(a) + W0 log(1-q_j(a)) - shrink a^2 (Newton steps, vectorised over replications)."""
    a = a0.copy(); Wt = W1 + W0
    for _ in range(iters):
        qa = _sig(lq[None, :] + a[:, None])
        g = (W1 - Wt * qa).sum(1) - 2 * shrink * a
        h = (Wt * qa * (1 - qa)).sum(1) + 2 * shrink
        a = np.clip(a + np.clip(g / h, -2, 2), -8, 8)
    return a


def sketch_marginal(jb, P, period_len, J):
    """Sketch of all units: counts (K,J) and sums of 1/p (R,K,J) per (period, bin)."""
    R, n = P.shape
    K = int(np.ceil(n / period_len))
    cnt = np.zeros((K, J)); sinv = np.zeros((R, K, J))
    for k in range(K):
        sl = slice(k * period_len, min(n, (k + 1) * period_len))
        oh = np.zeros((sl.stop - sl.start, J)); oh[np.arange(sl.stop - sl.start), jb[sl]] = 1.0
        cnt[k] = oh.sum(0); sinv[:, k] = (1.0 / P[:, sl]) @ oh
    return cnt, sinv


def sketch_sums(count_kj, sinv_kj, f_kj, v_kj, n):
    """Shore-side sums from the sketch: sum_t f_t and V_plug (independent draws: variance term (1-p)/p)."""
    sf = (f_kj * count_kj).sum((-1, -2))
    vp = (v_kj * (sinv_kj - count_kj)).sum((-1, -2)) / n ** 2
    return sf, vp


def prequential(y, S, P, jb, cal, period_len, shrink=1.0, mode='intercept', lam_floor=0.0):
    """Prequential estimator for fixed periods of period_len units (independent draws). Returns estimate, V-hat, V_plug and the intercepts."""
    R, n = S.shape; J = len(cal['q']); K = int(np.ceil(n / period_len))
    q, v = cal['q'], cal['v']; lq = _logit(q); q1q = q * (1 - q)
    cnt, sinv = sketch_marginal(jb, P, period_len, J)          # what the node uploads (marginal form)
    W1 = np.zeros((R, J)); W0 = np.zeros((R, J)); A = np.zeros(R); B = np.zeros(R)
    par = np.zeros((R, K)); a = np.zeros(R); lam = np.ones(R)
    resid = np.zeros(R); vres = np.zeros(R)
    for k in range(K):
        sl = slice(k * period_len, min(n, (k + 1) * period_len)); m = sl.stop - sl.start
        if mode == 'intercept':
            if k > 0: a = fit_intercept(lq, W1, W0, a, shrink)
            par[:, k] = a
            f = _sig(lq[jb[sl]][None, :] + a[:, None])
        elif mode == 'ppipp':
            if k > 0: lam = np.clip(np.where(B > 1e-300, A / np.maximum(B, 1e-300), 1.0), lam_floor, 1.0)
            par[:, k] = lam
            f = lam[:, None] * q[jb[sl]][None, :]
        else:
            f = np.broadcast_to(q[jb[sl]][None, :], (R, m))
        Sk, Pk, yk = S[:, sl], P[:, sl], y[sl][None, :]
        r = yk - f
        resid += (Sk * r / Pk).sum(1)
        vres += (Sk * (1 - Pk) * r ** 2 / Pk ** 2).sum(1)
        # update past-period statistics with period-k audited units (used only by periods > k)
        w = Sk / Pk
        oh = np.zeros((m, J)); oh[np.arange(m), jb[sl]] = 1.0
        if mode == 'intercept':
            W1 += (w * yk) @ oh; W0 += (w * (1 - yk)) @ oh
        elif mode == 'ppipp':
            fs = q[jb[sl]][None, :]; ww = Sk * (1 - Pk) / Pk ** 2
            A += (ww * yk * fs).sum(1); B += (ww * fs * fs).sum(1)
    # shore sums from the sketch only
    if mode == 'intercept':
        qa = _sig(lq[None, None, :] + par[:, :, None])                          # (R,K,J)
        fkj = qa; vkj = v[None, None, :] * qa * (1 - qa) / q1q[None, None, :]
    elif mode == 'ppipp':
        fkj = par[:, :, None] * q[None, None, :]
        vkj = v[None, None, :] + (1 - par[:, :, None]) ** 2 * q[None, None, :] ** 2
    else:
        fkj = np.broadcast_to(q[None, None, :], (R, K, J)); vkj = np.broadcast_to(v[None, None, :], (R, K, J))
    sf, vp = sketch_sums(cnt[None], sinv, fkj, vkj, n)
    return (sf + resid) / n, vres / n ** 2, vp, par


def preq_blocks(y, S, P, jb, cal, PID, shrink=0.3, vhat_fpc=False):
    """Prequential estimator with block-aligned periods PID (R,n) for EXCL and Madow; V-hat = sum xi r^2/p^2, V_plug = sum_kj v_kj sinv_kj / N^2. Returns estimate, V-hat and V_plug (R,)."""
    R, n = S.shape; J = len(cal['q']); K = int(PID.max()) + 1
    q, v = cal['q'], cal['v']; lq = _logit(q); q1q = q * (1 - q)
    idx = ((np.arange(R)[:, None] * K + PID.astype(np.int64)) * J + jb[None, :]).ravel(); sz = R * K * J
    w = (S / P)
    W1 = np.bincount(idx, (w * y[None, :]).ravel(), sz).reshape(R, K, J)
    W0 = np.bincount(idx, (w * (1 - y)[None, :]).ravel(), sz).reshape(R, K, J)
    cnt = np.bincount(idx, minlength=sz).reshape(R, K, J).astype(float)
    sinv = np.bincount(idx, (1.0 / P).ravel(), sz).reshape(R, K, J)
    a = np.zeros((R, K)); C1 = np.zeros((R, J)); C0 = np.zeros((R, J))
    for k in range(1, K):
        C1 += W1[:, k - 1]; C0 += W0[:, k - 1]
        a[:, k] = fit_intercept(lq, C1, C0, a[:, k - 1], shrink)
    au = np.take_along_axis(a, PID.astype(np.int64), 1)
    f = _sig(lq[jb][None, :] + au); r = y[None, :] - f
    resid = (S * r / P).sum(1)
    vhat = (S * r ** 2 / P ** 2).sum(1) if not vhat_fpc else (S * (1 - P) * r ** 2 / P ** 2).sum(1)
    qa = _sig(lq[None, None, :] + a[:, :, None]); vkj = v[None, None, :] * qa * (1 - qa) / q1q[None, None, :]
    sf = (qa * cnt).sum((1, 2))
    vp = (vkj * (sinv - (cnt if vhat_fpc else 0))).sum((1, 2)) / n ** 2
    return (sf + resid) / n, vhat / n ** 2, vp


def ht_est(y, S, P):
    n = len(y)
    est = (S * y / P).sum(1) / n
    var = (S * (1 - P) * y ** 2 / P ** 2).sum(1) / n ** 2
    return est, var


def ppi_est(f, y, S, P):
    n = len(y)
    r = y - f
    est = f.mean() + (S * r / P).sum(1) / n
    var = (S * (1 - P) * r ** 2 / P ** 2).sum(1) / n ** 2
    return est, var


def vplug_fixed(v, P):
    n = P.shape[1]
    return (v[None, :] * (1 - P) / P).sum(1) / n ** 2


def cm_width(theta, est, se):
    """Matched width: scale the half-widths so that the empirical coverage is the smallest value >= 0.95 over the R runs; returns the mean full width after scaling."""
    R = len(est); se = np.asarray(se); err = np.abs(est - theta)
    z = np.where(se > 0, err / np.where(se > 0, se, 1), np.where(err <= 1e-15, 0.0, np.inf))
    k = int(np.ceil(0.95 * R - 1e-9)); zk = np.sort(z)[k - 1]
    if not np.isfinite(zk): return None
    return float(2 * zk * se.mean())


def block_v1v2(xr, NEW, rec_unit, n):
    """Block-total variance estimates. xr = xi r / p (R,n); NEW block starts; rec_unit recording index of every unit. V1 = N^-2 sum_k T_k^2 (equals V-hat). V2 = N^-2 (1/2) sum_k w_k with w_k = (T_k - T_{k+1})^2 for a block with a successor in its recording, w_k = (T_{k-1} - T_k)^2 for the last block of a recording with at least two blocks, w_k = 2 T_k^2 for a single-block recording (successive-difference estimate). Returns V1, V2 and the number of blocks."""
    R = xr.shape[0]; bid = np.cumsum(NEW, 1) - 1; nb = bid[:, -1] + 1; Kb = int(nb.max())
    Z = np.bincount((np.arange(R)[:, None] * Kb + bid).ravel(), xr.ravel(), R * Kb).reshape(R, Kb)
    rows, cols = np.nonzero(NEW)
    recb = np.full((R, Kb), -1, np.int64); recb[rows, bid[rows, cols]] = rec_unit[cols]
    valid = recb >= 0
    nxt = np.zeros((R, Kb), bool); nxt[:, :-1] = recb[:, :-1] == recb[:, 1:]
    prv = np.zeros((R, Kb), bool); prv[:, 1:] = recb[:, 1:] == recb[:, :-1]
    Zn = np.zeros_like(Z); Zn[:, :-1] = Z[:, 1:]; Zp = np.zeros_like(Z); Zp[:, 1:] = Z[:, :-1]
    w = np.where(nxt, (Z - Zn) ** 2, np.where(prv, (Zp - Z) ** 2, 2.0 * Z ** 2)) * valid
    return (Z ** 2 * valid).sum(1) / n ** 2, w.sum(1) / 2.0 / n ** 2, nb


def coverage(theta, est, se):
    """Fraction of replications whose interval est +- Z95 se contains theta."""
    return float(((est - Z95 * se <= theta) & (theta <= est + Z95 * se)).mean())
