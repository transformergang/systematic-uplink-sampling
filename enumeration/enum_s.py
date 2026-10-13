#!/usr/bin/env python
"""Exact enumeration / linear-programming verification of Proposition 1 and Corollary 1 (checks E1-E7).

Conventions
- Unit t (0-based) is bit t of the state index x in [0, 2^N): outcome vector xi_t = (x >> t) & 1.
- A partition is a list of half-open unit ranges [(lo, hi), ...]; the history of the block starting at lo is the
  low `lo` bits of x (x & (2^lo - 1)).
- p_hist[t] is an array over the histories of t's block start (length 2^lo): p_t(h).  Open loop = constant array.
- Run:  python enum_s.py [--n-inst 200] [--seed S] [--only E1,E3] [--smoke]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, 'results')
PLUGIN_DIR = os.path.join(os.path.dirname(HERE), 'code', 'python')

LP_TOL = 1e-7        # tolerance for quantities that come out of an LP solver
EXACT_TOL = 1e-10    # tolerance for closed-form enumeration
MASS_TOL = 1e-12     # slack when testing "block mass <= 1"
MAX_DUMP = 8         # max failing instances dumped per check


# --------------------------------------------------------------------------- bookkeeping
class Tracker:
    """Collects per-item max errors / flags, and the failing instances (dumped to JSON)."""

    def __init__(self):
        self.items = {}
        self.failing = []
        self.notes = {}

    def err(self, name, err, tol, inst=None):
        """Numeric item: fails if err > tol (err should be >= 0; for inequalities pass the violation)."""
        err = float(err)
        if not np.isfinite(err):
            err = float('inf')
        d = self.items.setdefault(name, dict(kind='err', n=0, n_fail=0, max_err=0.0, tol=tol))
        d['n'] += 1
        d['max_err'] = max(d['max_err'], err)
        if err > tol:
            d['n_fail'] += 1
            self._dump(name, inst, err=err, tol=tol)

    def flag(self, name, ok, inst=None, strict=True):
        """Boolean item.  strict=False: counted and exemplified in `notes`, but never makes the check fail."""
        d = self.items.setdefault(name, dict(kind='flag' if strict else 'info_flag', n=0, n_fail=0))
        d['n'] += 1
        if not ok:
            d['n_fail'] += 1
            if strict:
                self._dump(name, inst)
            else:
                ex = self.notes.setdefault(name + '_examples', [])
                if len(ex) < 5:
                    ex.append(inst)

    def _dump(self, name, inst, **kw):
        if len(self.failing) < MAX_DUMP:
            self.failing.append(dict(item=name, instance=inst, **kw))

    def result(self, n_instances, **info):
        errs = [d['max_err'] for d in self.items.values() if d['kind'] == 'err']
        strict_fail = sum(d['n_fail'] for d in self.items.values() if d['kind'] in ('err', 'flag'))
        return dict(n_instances=int(n_instances), max_abs_err=max(errs) if errs else 0.0,
                    pass_=bool(strict_fail == 0), n_strict_failures=int(strict_fail), failing=self.failing,
                    items=self.items, notes=self.notes, info=info)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k).rstrip('_') if k == 'pass_' else str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


# --------------------------------------------------------------------------- bit / partition / event helpers
_BITS = {}


def get_bits(N):
    """(2^N, N) int array, column t = outcome of unit t in each state."""
    if N not in _BITS:
        _BITS[N] = ((np.arange(1 << N)[:, None] >> np.arange(N)[None, :]) & 1).astype(np.int8)
    return _BITS[N]


def excl_partition(p, forced):
    """EXCL block rule (same float arithmetic as samplers.simulate_excl_openloop): new block when A_{t-1}+p_t > 1 or forced."""
    blocks, start, A = [], 0, 0.0
    for t in range(len(p)):
        if t > 0 and (forced[t] or A + p[t] > 1.0):
            blocks.append((start, t))
            start, A = t, 0.0
        A += p[t]
    blocks.append((start, len(p)))
    return blocks


def event_parts(blocks, e):
    """[(block index, units of e in that block)] for the blocks that meet e."""
    es = set(e)
    out = []
    for k, (lo, hi) in enumerate(blocks):
        units = [t for t in range(lo, hi) if t in es]
        if units:
            out.append((k, units))
    return out


def union_bound_Ue(blocks, p, e):
    """U_e = 1 - prod_k (1 - min(1, m_k)) and the list of part masses m_k."""
    ms = [float(sum(p[t] for t in units)) for _, units in event_parts(blocks, e)]
    return 1.0 - float(np.prod([1.0 - min(1.0, m) for m in ms])), ms


def event_mask(e):
    return sum(1 << int(t) for t in e)


def random_event(rng, N, p, kind):
    """kind: 'contig' | 'subset' | 'heavy' (random order prefix with sum p > 1, if possible)."""
    if kind == 'contig':
        lo = int(rng.integers(N))
        hi = int(rng.integers(lo, N))
        return list(range(lo, hi + 1))
    if kind == 'subset':
        while True:
            m = rng.random(N) < 0.5
            if m.any():
                return [int(t) for t in np.flatnonzero(m)]
    if kind == 'heavy':
        order, s, e = rng.permutation(N), 0.0, []
        for t in order:
            e.append(int(t))
            s += p[t]
            if s > 1.0 and len(e) >= 2:
                break
        return sorted(e)
    raise ValueError(kind)


def block_subset_event(rng, blocks):
    """Random subset (>= 2 units) inside one block of size >= 2, or None."""
    big = [b for b in blocks if b[1] - b[0] >= 2]
    if not big:
        return None
    lo, hi = big[int(rng.integers(len(big)))]
    k = int(rng.integers(2, hi - lo + 1))
    return sorted(int(t) for t in rng.choice(np.arange(lo, hi), size=k, replace=False))


# --------------------------------------------------------------------------- instance generation
def gen_instance(rng, N, mode, force_prob=0.3):
    """Open-loop instance (p, blocks).  mode: 'std' (block masses <= 1 by the EXCL rule), 'mass1' (one block of mass
    exactly 1), 'big' (one block with mass in (1, 1.6]).  Falls back to 'std' if no suitable block exists."""
    p = rng.uniform(0.02, 0.98, N)
    forced = rng.random(N) < force_prob
    blocks = excl_partition(p, forced)
    if mode == 'std':
        return dict(N=N, p=p, blocks=blocks, mode='std')
    if all(hi - lo < 2 for lo, hi in blocks) and len(blocks) >= 2:       # merge two blocks so a size>=2 block exists
        blocks = [(blocks[0][0], blocks[1][1])] + blocks[2:]
    cand = [b for b in blocks if 2 <= b[1] - b[0] <= (7 if mode == 'mass1' else 3)]
    if not cand:
        return dict(N=N, p=p, blocks=blocks, mode='std')
    lo, hi = cand[int(rng.integers(len(cand)))]
    s = hi - lo
    if mode == 'mass1':
        p[lo:hi] = 0.02 + (1.0 - 0.02 * s) * rng.dirichlet(np.ones(s))   # sums to 1, each in [0.02, 0.98]
        return dict(N=N, p=p, blocks=blocks, mode='mass1')
    for _ in range(500):                                                  # 'big'
        q = rng.uniform(1.05, 1.6) * rng.dirichlet(np.ones(s))
        if q.min() >= 0.02 and q.max() <= 0.98:
            p[lo:hi] = q
            return dict(N=N, p=p, blocks=blocks, mode='big')
    return dict(N=N, p=p, blocks=blocks, mode='std')


def instance_dump(inst, **extra):
    d = dict(N=inst['N'], p=inst['p'], blocks=inst['blocks'], mode=inst.get('mode'))
    d.update(extra)
    return d


def block_masses_ok(blocks, p_hist):
    """True iff every block has mass <= 1 (up to MASS_TOL) for every history."""
    for lo, hi in blocks:
        mass = np.sum([p_hist[t] for t in range(lo, hi)], axis=0)
        if mass.max() > 1.0 + MASS_TOL:
            return False
    return True


# --------------------------------------------------------------------------- p_t(h) tables, constraints, laws
def hist_bits(lo, hide_last):
    """Number of low bits visible at a block start; 'delayed information' hides the last unit before the block."""
    return lo - 1 if (hide_last and lo >= 1) else lo


def open_p_hist(N, blocks, p, hide_last=False):
    ph = [None] * N
    for lo, hi in blocks:
        for t in range(lo, hi):
            ph[t] = np.full(1 << hist_bits(lo, hide_last), float(p[t]))
    return ph


def closed_p_hist(rng, N, blocks):
    """Random p_t(h) in [0.02, 0.98], rescaled where needed so that every block mass is <= 1 for every history."""
    ph = [None] * N
    for lo, hi in blocks:
        arr = rng.uniform(0.02, 0.98, (hi - lo, 1 << lo))
        mass = arr.sum(0)
        arr = arr * np.where(mass > 1.0, rng.uniform(0.6, 1.0, mass.shape) / mass, 1.0)
        for j, t in enumerate(range(lo, hi)):
            ph[t] = arr[j]
    return ph


def conditional_constraints(N, blocks, p_hist, hide_last=False):
    """A_eq, b_eq of Prop. 1 (conditional marginals) plus sum(P)=1:  for every block, history h, unit t in block:
         sum_{x: hist(x)=h, x_t=1} P(x) - p_t(h) * sum_{x: hist(x)=h} P(x) = 0."""
    idx = np.arange(1 << N)
    bits = get_bits(N)
    rows, cols, vals, nrow = [], [], [], 0
    for lo, hi in blocks:
        keep = hist_bits(lo, hide_last)
        hist = idx & ((1 << keep) - 1)
        for t in range(lo, hi):
            assert len(p_hist[t]) == 1 << keep
            rows.append(nrow + hist)
            cols.append(idx)
            vals.append(bits[:, t] - p_hist[t][hist])
            nrow += 1 << keep
    rows.append(np.full(1 << N, nrow))
    cols.append(idx)
    vals.append(np.ones(1 << N))
    A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(nrow + 1, 1 << N)).tocsr()
    b = np.zeros(nrow + 1)
    b[-1] = 1.0
    return A, b


def marginal_constraints(N, p):
    """Only the unconditional marginals P(xi_t = 1) = p_t, and sum(P) = 1."""
    A = sp.csr_matrix(np.vstack([get_bits(N).T.astype(float), np.ones((1, 1 << N))]))
    return A, np.concatenate([p, [1.0]])


def solve_lp(c, A, b):
    res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    if res.status != 0:
        raise RuntimeError('LP failed: ' + str(res.message))
    return np.clip(res.x, 0.0, None), float(res.fun)


def lp_max_hit(N, A, b, e):
    """max P(e hit) over P with A P = b, P >= 0.  Returns (optimum, argmax P)."""
    hit = ((np.arange(1 << N) & event_mask(e)) != 0).astype(float)
    P, fun = solve_lp(-hit, A, b)
    return -fun, P


def excl_law(N, blocks, p_hist):
    """Exact EXCL law (vector over 2^N): in block k, given the realised history h, either one unit t uploads
    (prob p_t(h)) or none (prob 1 - m_k(h)); blocks are conditionally independent given the past."""
    w = np.ones(1)
    for lo, hi in blocks:
        ph = np.array([p_hist[t] for t in range(lo, hi)])             # (s, 2^lo)
        if ph.sum(0).max() > 1.0 + MASS_TOL:
            raise ValueError('EXCL needs block mass <= 1')
        s = hi - lo
        new = np.zeros(1 << hi)
        h = np.arange(1 << lo)
        new[h] = w * (1.0 - ph.sum(0))                                # no upload in the block
        for j in range(s):
            new[h + (1 << (lo + j))] = w * ph[j]
        w = new
    return w


def indep_law(p):
    P = np.ones(1)
    for pt in p:
        P = np.concatenate([P * (1.0 - pt), P * pt])
    return P


def madow_law(p):
    """Exact law of single-start systematic sampling: u ~ U(0,1], unit t selected iff some u+j (j>=0 integer) lies in
    (C_{t-1}, C_t], C = global cumulative mass.  (0,1] is cut at the fractional parts of C; the outcome is constant on
    each piece, so it is evaluated at the piece midpoint."""
    N = len(p)
    C = np.concatenate([[0.0], np.cumsum(p)])
    cuts = np.unique(np.concatenate([[0.0, 1.0], np.mod(C, 1.0)]))
    law = np.zeros(1 << N)
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a <= 0:
            continue
        u = 0.5 * (a + b)
        cnt = np.floor(C[1:] - u) - np.floor(C[:-1] - u)              # number of points u+j in (C_{t-1}, C_t]
        assert cnt.max() <= 1
        law[int(sum(1 << t for t in range(N) if cnt[t] > 0))] += b - a
    return law


def decompose(P, N, blocks, p_hist, e):
    """E[sum_k M_k], P(e hit), R_e, W_e of Prop. 1 (decomposition) for the law P (vector over 2^N)."""
    idx = np.arange(1 << N)
    bits = get_bits(N)
    D = np.ones(1 << N)
    sumM = np.zeros(1 << N)
    Rterm = np.zeros(1 << N)
    Wterm = np.zeros(1 << N)
    for k, units in event_parts(blocks, e):
        lo = blocks[k][0]
        hist = idx & ((1 << lo) - 1)
        M = sum(p_hist[t][hist] for t in units)
        Nk = sum(bits[:, t].astype(int) for t in units)
        sumM += M
        Rterm += (1 - D) * M
        Wterm += D * np.maximum(Nk - 1, 0)
        D = D * (Nk == 0)
    hit = ((idx & event_mask(e)) != 0).astype(float)
    return dict(EM=float(P @ sumM), hit=float(P @ hit), R=float(P @ Rterm), W=float(P @ Wterm))


def cond_violation(P, N, lo, t):
    """Per history h of the block start `lo`: (P(h), P(xi_t=1 | h)) where P(h) > 0."""
    idx = np.arange(1 << N)
    hist = idx & ((1 << lo) - 1)
    Ph = np.bincount(hist, weights=P, minlength=1 << lo)
    Pt = np.bincount(hist, weights=P * get_bits(N)[:, t], minlength=1 << lo)
    ok = Ph > 1e-14
    q = np.zeros_like(Ph)
    q[ok] = Pt[ok] / Ph[ok]
    return Ph, q, ok


# --------------------------------------------------------------------------- E1
def check_E1(n_inst, rng):
    """LP max P(e hit) under conditional exactness (open loop) = U_e; EXCL attains it with W_e = 0."""
    tr = Tracker()
    n_ev, n_heavy, n_big_attained, n_big = 0, 0, 0, 0
    for i in range(n_inst):
        mode = ['std', 'std', 'std', 'mass1', 'big'][i % 5]
        inst = gen_instance(rng, int(rng.integers(3, 8)), mode)
        N, p, blocks = inst['N'], inst['p'], inst['blocks']
        ph = open_p_hist(N, blocks, p)
        A, b = conditional_constraints(N, blocks, ph)
        masses_ok = block_masses_ok(blocks, ph)
        if masses_ok:
            law = excl_law(N, blocks, ph)
            tr.err('excl_law_feasible', np.abs(A @ law - b).max(), EXACT_TOL, instance_dump(inst))
        for kind in ('contig', 'subset', 'heavy'):
            e = random_event(rng, N, p, kind)
            U, ms = union_bound_Ue(blocks, p, e)
            opt, _ = lp_max_hit(N, A, b, e)
            dump = instance_dump(inst, e=e, lp=opt, U=U, kind=kind)
            n_ev += 1
            n_heavy += float(sum(p[t] for t in e)) > 1.0
            if masses_ok:
                tr.err('lp_minus_U_abs', abs(opt - U), LP_TOL, dump)           # bound AND attainment
                st = decompose(law, N, blocks, ph, e)
                tr.err('excl_hit_minus_U', abs(st['hit'] - U), EXACT_TOL, dump)
                tr.err('excl_W_abs', abs(st['W']), EXACT_TOL, dump)
                if len(ms) == 2 and sum(ms) <= 1.0:                            # loss = m1*m2 for a two-block split
                    tr.err('loss_equals_m1m2', abs(min(1.0, sum(ms)) - U - ms[0] * ms[1]), EXACT_TOL, dump)
            else:                                                              # inequality only; attainment not claimed
                tr.err('big_lp_exceeds_U', max(0.0, opt - U), LP_TOL, dump)
                n_big += 1
                n_big_attained += abs(opt - U) <= LP_TOL
    return tr.result(n_inst, n_events=n_ev, n_events_sum_p_gt_1=int(n_heavy), n_events_on_big_block_instances=n_big,
                     n_big_events_where_bound_attained_anyway=int(n_big_attained))


# --------------------------------------------------------------------------- E2
def check_E2(n_inst, rng):
    """E[sum M] - P(hit) = R_e + W_e (and R, W >= 0) on random LP vertices, open and closed loop."""
    tr = Tracker()
    n_nontrivial = 0
    for i in range(n_inst):
        closed = bool(i % 2)
        inst = gen_instance(rng, int(rng.integers(3, 8)), ['std', 'mass1', 'big'][(i // 2) % 3])
        N, p, blocks = inst['N'], inst['p'], inst['blocks']
        ph = closed_p_hist(rng, N, blocks) if closed else open_p_hist(N, blocks, p)
        A, b = conditional_constraints(N, blocks, ph)
        P, _ = solve_lp(rng.normal(size=1 << N), A, b)                        # random vertex
        feas = np.abs(A @ P - b).max()
        for kind in ('contig', 'subset', 'heavy'):
            e = random_event(rng, N, p, kind)
            st = decompose(P, N, blocks, ph, e)
            dump = instance_dump(inst, e=e, closed=closed, stats=st, feas=feas)
            tr.err('feasible', feas, LP_TOL, dump)
            tr.err('identity_abs', abs(st['EM'] - st['hit'] - st['R'] - st['W']), LP_TOL, dump)
            tr.err('R_nonneg_violation', max(0.0, -st['R']), LP_TOL, dump)
            tr.err('W_nonneg_violation', max(0.0, -st['W']), LP_TOL, dump)
            n_nontrivial += (st['R'] > 1e-6 and st['W'] > 1e-6)
    return tr.result(n_inst, n_events=3 * n_inst, n_events_with_R_and_W_both_positive=int(n_nontrivial))


# --------------------------------------------------------------------------- E3 / E4 (shared design generation)
def crossing_event(rng, blocks, N):
    """Contiguous event containing the last unit of some block and the first unit of the next one."""
    s = blocks[int(rng.integers(1, len(blocks)))][0]
    return list(range(int(rng.integers(0, s)), int(rng.integers(s, N)) + 1))


def gen_crossing_instance(rng):
    N = int(rng.integers(3, 8))
    p = rng.uniform(0.02, 0.98, N)
    forced = rng.random(N) < 0.3
    blocks = excl_partition(p, forced)
    if len(blocks) < 2:
        forced[N // 2] = True
        blocks = excl_partition(p, forced)
    return dict(N=N, p=p, blocks=blocks, mode='std')


def e3_designs(n_inst, rng):
    """Instances + marginal-only LP optimiser for a boundary-crossing contiguous event (shared by E3(i) and E4)."""
    for _ in range(n_inst):
        inst = gen_crossing_instance(rng)
        e = crossing_event(rng, inst['blocks'], inst['N'])
        A, b = marginal_constraints(inst['N'], inst['p'])
        opt, P = lp_max_hit(inst['N'], A, b, e)
        yield inst, e, opt, P


def check_E3(n_inst, rng):
    tr = Tracker()
    # (i) marginals only
    n_strict = 0
    for inst, e, opt, _ in e3_designs(n_inst, rng):
        N, p, blocks = inst['N'], inst['p'], inst['blocks']
        target = min(1.0, float(sum(p[t] for t in e)))
        U, ms = union_bound_Ue(blocks, p, e)
        mad = madow_law(p)
        mad_hit = float(mad[(np.arange(1 << N) & event_mask(e)) != 0].sum())
        dump = instance_dump(inst, e=e, lp=opt, target=target, U=U, madow=mad_hit)
        tr.err('i_lp_minus_min1sum', abs(opt - target), LP_TOL, dump)
        tr.err('i_madow_minus_min1sum', abs(mad_hit - target), EXACT_TOL, dump)
        tr.err('i_madow_marginals', np.abs(get_bits(N).T @ mad - p).max(), EXACT_TOL, dump)
        tr.err('i_gap_negative_violation', max(0.0, U - opt), LP_TOL, dump)
        if len(ms) >= 2:                                                     # event meets >= 2 blocks: strict gap
            tr.flag('i_strict_gap_gt_0', opt - U > 1e-9, dump)
            n_strict += 1
    # (ii) delayed information: outcome of the last unit before each block start hidden
    best_gap, best_inst = -1.0, None
    for _ in range(n_inst):
        inst = gen_crossing_instance(rng)
        N, p, blocks = inst['N'], inst['p'], inst['blocks']
        e = crossing_event(rng, blocks, N)
        A, b = conditional_constraints(N, blocks, open_p_hist(N, blocks, p, hide_last=True), hide_last=True)
        opt, _ = lp_max_hit(N, A, b, e)
        U, _ = union_bound_Ue(blocks, p, e)
        target = min(1.0, float(sum(p[t] for t in e)))
        dump = instance_dump(inst, e=e, lp=opt, U=U, target=target)
        tr.err('ii_below_U_violation', max(0.0, U - opt), LP_TOL, dump)       # EXCL stays feasible
        tr.err('ii_above_unionbound_violation', max(0.0, opt - target), LP_TOL, dump)
        if opt - U > best_gap:
            best_gap, best_inst = opt - U, dump
    tr.flag('ii_exists_instance_with_gap_gt_1e-6', best_gap > 1e-6, best_inst)
    return tr.result(n_inst, n_i_events_meeting_ge2_blocks=n_strict, ii_max_gap_LP_minus_U=best_gap,
                     ii_best_instance=best_inst)


def theta_hat_vectors(N, p, y, f_cols):
    """theta_hat_f(x) for every state x; f_cols (2^N, N) are the predictor values in each state."""
    bits = get_bits(N)
    return (f_cols + bits * (y[None, :] - f_cols) / p[None, :]).sum(1) / N


def check_E4(n_inst, rng):
    """Designs from E3(i) violating conditional exactness are biased for the PROOF-idea predictor."""
    tr = Tracker()
    n_viol, n_nonviol = 0, 0
    for inst, e, opt, P in e3_designs(n_inst, rng):
        N, p, blocks = inst['N'], inst['p'], inst['blocks']
        idx = np.arange(1 << N)
        # worst violated (t, h): weighted by P(h)
        worst, arg = 0.0, None
        for lo, hi in blocks[1:]:
            for t in range(lo, hi):
                Ph, q, ok = cond_violation(P, N, lo, t)
                dev = np.where(ok, Ph * np.abs(q - p[t]), 0.0)
                if dev.max() > worst:
                    worst, arg = float(dev.max()), (lo, t, Ph, q, ok)
        if worst <= 1e-6:                                                    # (not expected: LP optimum exceeds U_e)
            n_nonviol += 1
            continue
        n_viol += 1
        lo, t, Ph, q, ok = arg
        theta = 1.0 / N                                                      # y_t = 1 at the single unit t
        y = np.zeros(N)
        y[t] = 1.0
        hist = idx & ((1 << lo) - 1)
        excl = excl_law(N, blocks, open_p_hist(N, blocks, p))                # control: valid design -> zero bias
        biases = {}
        for name, A_set in (('A_lt', ok & (q < p[t])), ('A_gt', ok & (q > p[t]))):
            f_cols = np.tile(y, (1 << N, 1))                                 # f_s = y_s for s != t
            f_cols[:, t] = A_set[hist].astype(float)                         # f_t = 1_A(h), h = history at t's block start
            th = theta_hat_vectors(N, p, y, f_cols)
            biases[name] = float(P @ th - theta)
            dump = instance_dump(inst, e=e, t=t, block_start=lo, worst_violation=worst, bias=biases[name], A=name)
            tr.flag('bias_abs_gt_1e-9_' + name, abs(biases[name]) > 1e-9, dump)
            tr.err('control_excl_bias_' + name, abs(excl @ th - theta), EXACT_TOL, dump)
        tr.items.setdefault('min_abs_bias', dict(kind='stat', value=np.inf))
        tr.items['min_abs_bias']['value'] = min(tr.items['min_abs_bias']['value'], min(abs(v) for v in biases.values()))
    tr.flag('at_least_one_violating_design', n_viol >= 1)
    return tr.result(n_inst, n_violating_designs=n_viol, n_designs_not_violating=n_nonviol)


# --------------------------------------------------------------------------- E5
def dprime_law(pa, pb, g, x):
    """Joint law of (xi_a, xi_b, xi_c) (bits 0,1,2): P(xi_a=xi_b=1)=x, then xi_c ~ Bernoulli(g[#uploads in B_1])."""
    P = np.zeros(8)
    first = {(1, 1): x, (1, 0): pa - x, (0, 1): pb - x, (0, 0): 1.0 - pa - pb + x}
    for (xa, xb), w in first.items():
        s = xa | (xb << 1)
        P[s] += w * (1.0 - g[xa + xb])
        P[s | 4] += w * g[xa + xb]
    return P


def e5_case(pa, pb, g, x, tr, tag, g2_equals_g1):
    blocks = [(0, 2), (2, 3)]
    ph = [np.array([pa]), np.array([pb]), np.array([g[0], g[1], g[1], g[2]])]  # p_c(h), h = (a, b) bits
    e = [0, 2]                                                               # e = {a, c}
    Pe, Pd = excl_law(3, blocks, ph), dprime_law(pa, pb, g, x)
    se, sd = decompose(Pe, 3, blocks, ph, e), decompose(Pd, 3, blocks, ph, e)
    A, b = conditional_constraints(3, blocks, ph)
    dump = dict(tag=tag, pa=pa, pb=pb, g=g, x=x, P_excl=se['hit'], P_dprime=sd['hit'])
    tr.err('dprime_feasible_residual', np.abs(A @ Pd - b).max(), EXACT_TOL, dump)
    tr.err('excl_feasible_residual', np.abs(A @ Pe - b).max(), EXACT_TOL, dump)
    tr.err('dprime_nonneg_violation', max(0.0, -Pd.min()), EXACT_TOL, dump)
    tr.err('hit_diff_vs_x(g0-g1)', abs(sd['hit'] - se['hit'] - x * (g[0] - g[1])), EXACT_TOL, dump)
    tr.err('Emass_diff_vs_general_x(g0-2g1+g2)', abs(sd['EM'] - se['EM'] - x * (g[0] - 2 * g[1] + g[2])), EXACT_TOL, dump)
    if g2_equals_g1:
        tr.err('Emass_diff_vs_claim_x(g0-g1)', abs(sd['EM'] - se['EM'] - x * (g[0] - g[1])), EXACT_TOL, dump)
    return se, sd


def check_E5(n_inst, rng):
    tr = Tracker()
    se, sd = e5_case(0.3, 0.4, [0.6, 0.2, 0.2], 0.25, tr, 'fixed Prop. 1 (closed-loop counterexample) numbers (g(2):=g(1))', True)
    tr.err('fixed_P_excl_is_0.56', abs(se['hit'] - 0.56), EXACT_TOL)
    tr.err('fixed_P_dprime_is_0.66', abs(sd['hit'] - 0.66), EXACT_TOL)
    for i in range(n_inst):
        while True:
            pa, pb = rng.uniform(0.02, 0.9, 2)
            if pa + pb <= 1.0:
                break
        g0 = rng.uniform(0.2, 0.98)
        g1 = rng.uniform(0.02, g0)
        x = rng.uniform(0.01, 1.0) * min(pa, pb)
        e5_case(pa, pb, [g0, g1, g1], x, tr, 'grid g2=g1', True)
        e5_case(pa, pb, [g0, g1, rng.uniform(0.02, 0.98)], x, tr, 'grid g2 random', False)
    return tr.result(1 + n_inst, fixed=dict(P_excl=se['hit'], P_dprime=sd['hit'], diff=sd['hit'] - se['hit']))


# --------------------------------------------------------------------------- E6
def elem_sym(v):
    """sigma_1, sigma_2, sigma_3 of the values v."""
    s = [1.0, 0.0, 0.0, 0.0]
    for x in v:
        for j in (3, 2, 1):
            s[j] += s[j - 1] * x
    return s[1], s[2], s[3]


def variance(P, vals):
    m = P @ vals
    return float(P @ (vals - m) ** 2)


def e6_events(rng, N, p, blocks):
    evs = [random_event(rng, N, p, k) for k in ('contig', 'subset', 'heavy', 'subset')]
    w = block_subset_event(rng, blocks)
    if w is not None:
        evs.append(w)
    ones = [t for t in range(N) if p[t] >= 1.0]
    if ones:
        evs.append(sorted(set(ones[:1]) | set(random_event(rng, N, p, 'subset'))))
        multi = [b for b in blocks if b[1] - b[0] >= 2 and not (b[0] <= ones[0] < b[1])]
        if multi:
            lo, hi = multi[int(rng.integers(len(multi)))]
            evs.append(sorted({ones[0]} | set(range(lo, hi))))               # a p=1 unit plus a multi-unit part
    return evs


def e6_laws_checks(tr, inst):
    """Parts of Cor. 1 that concern the whole design: properties of the whole design (no event)."""
    N, p, blocks, rng = inst['N'], inst['p'], inst['blocks'], inst['rng']
    bits = get_bits(N)
    idx = np.arange(1 << N)
    ph = open_p_hist(N, blocks, p)
    PE, PI = excl_law(N, blocks, ph), indep_law(p)
    dump = instance_dump(inst)
    ms = np.array([p[lo:hi].sum() for lo, hi in blocks])
    # (a)
    tr.err('a_excl_sum_to_1', abs(PE.sum() - 1), EXACT_TOL, dump)
    tr.err('a_marginals_excl', np.abs(bits.T @ PE - p).max(), EXACT_TOL, dump)
    tr.err('a_marginals_indep', np.abs(bits.T @ PI - p).max(), EXACT_TOL, dump)
    A, b = conditional_constraints(N, blocks, ph)
    tr.err('a_conditional_exactness_excl', np.abs(A @ PE - b).max(), EXACT_TOL, dump)
    S_all = bits.sum(1).astype(float)
    tr.err('a_E_sum_xi', max(abs(PE @ S_all - p.sum()), abs(PI @ S_all - p.sum())), EXACT_TOL, dump)
    vE, vI = variance(PE, S_all), variance(PI, S_all)
    pair = sum(p[lo:hi].sum() ** 2 - (p[lo:hi] ** 2).sum() for lo, hi in blocks)   # sum_k sum_{s != t} p_s p_t
    tr.err('a_VarE_formula', abs(vE - (ms * (1 - ms)).sum()), EXACT_TOL, dump)
    tr.err('a_VarI_formula', abs(vI - (p * (1 - p)).sum()), EXACT_TOL, dump)
    tr.err('a_VarE_eq_VarI_minus_pairs', abs(vE - (vI - pair)), EXACT_TOL, dump)
    # (b) random history-dependent predictor f_t(h) (values in [0,1]), random labels
    y = (rng.random(N) < 0.5).astype(float)
    f_cols = np.zeros((1 << N, N))
    for lo, hi in blocks:
        for t in range(lo, hi):
            f_cols[:, t] = rng.uniform(0, 1, 1 << lo)[idx & ((1 << lo) - 1)]
    r = y[None, :] - f_cols
    th = theta_hat_vectors(N, p, y, f_cols)
    theta = y.mean()
    tr.err('b_unbiased', abs(PE @ th - theta), EXACT_TOL, dump)
    Vhat = (bits * r ** 2 / p[None, :] ** 2).sum(1) / N ** 2
    blocksum2 = sum(r[:, lo:hi].sum(1) ** 2 for lo, hi in blocks) / N ** 2
    lhs = PE @ Vhat - variance(PE, th)
    tr.err('b_EVhat_minus_Var', abs(lhs - PE @ blocksum2), EXACT_TOL, dump)
    tr.err('b_nonneg_violation', max(0.0, -lhs), EXACT_TOL, dump)
    # (d) deterministic predictor: residual r in [-1, 1]
    rd = rng.uniform(-1, 1, N)
    C = sum(rd[lo:hi].sum() ** 2 - (rd[lo:hi] ** 2).sum() for lo, hi in blocks)
    S = float((rd ** 2 * (1 - p) / p).sum())
    dev = (bits * (rd / p)[None, :]).sum(1) / N                              # theta_hat - theta - const; variance only
    VE, VI = variance(PE, dev), variance(PI, dev)
    tr.err('d_VarE_minus_VarI_eq_minusC', abs(VE - VI + C / N ** 2), EXACT_TOL, dump)
    tr.err('d_VarI_eq_S', abs(VI - S / N ** 2), EXACT_TOL, dump)
    if S > 1e-9:
        tr.err('d_ratio_eq_1_minus_C_over_S', abs(VE / VI - (1 - C / S)), EXACT_TOL, dump)
    VhE = (bits * (rd ** 2 / p ** 2)[None, :]).sum(1) / N ** 2
    VhI = (bits * ((1 - p) * rd ** 2 / p ** 2)[None, :]).sum(1) / N ** 2
    tr.err('d_EVhatI_unbiased_under_I', abs(PI @ VhI - VI), EXACT_TOL, dump)
    tr.err('d_EEVhatE_minus_EIVhatI', abs(PE @ VhE - PI @ VhI - (rd ** 2).sum() / N ** 2), EXACT_TOL, dump)
    return PE, PI


def e6_event_checks(tr, inst, PE, PI, e, kind):
    """Part of Cor. 1 that concerns one event for one event."""
    N, p, blocks = inst['N'], inst['p'], inst['blocks']
    hit = ((np.arange(1 << N) & event_mask(e)) != 0)
    parts = event_parts(blocks, e)
    m = np.array([sum(p[t] for t in u) for _, u in parts])
    a = np.array([np.prod([1 - p[t] for t in u]) for _, u in parts])
    g = a - (1 - m)
    dump = instance_dump(inst, e=e, kind=kind)
    PEh, PIh = float(PE[hit].sum()), float(PI[hit].sum())
    G = PEh - PIh
    tr.err('c_PE_formula', abs(PEh - (1 - np.prod(1 - m))), EXACT_TOL, dump)
    tr.err('c_PI_formula', abs(PIh - (1 - np.prod([1 - p[t] for t in e]))), EXACT_TOL, dump)
    tr.err('c_dominance_violation', max(0.0, PIh - PEh), EXACT_TOL, dump)
    lower = sum(g[k] * np.prod(np.delete(1 - m, k)) for k in range(len(parts)))
    tr.err('c_bracket_lower_violation', max(0.0, lower - G), EXACT_TOL, dump)
    tr.err('c_bracket_upper_violation', max(0.0, G - g.sum()), EXACT_TOL, dump)
    for (k, u), gk in zip(parts, g):
        s1, s2, s3 = elem_sym([p[t] for t in u])
        tr.err('c_gk_lower_violation', max(0.0, (s2 - s3) - gk), EXACT_TOL, dump)
        tr.err('c_gk_upper_violation', max(0.0, gk - s2), EXACT_TOL, dump)
    if len(parts) == 1:
        s1, s2, s3 = elem_sym([p[t] for t in e])
        tr.err('c_within_block_bracket_violation', max(0.0, (s2 - s3) - G, G - s2), EXACT_TOL, dump)
    sum_p = float(sum(p[t] for t in e))
    tr.err('c_union_gap_formula', abs((min(1.0, sum_p) - PEh) - (min(1.0, m.sum()) - 1 + np.prod(1 - m))),
           EXACT_TOL, dump)
    tr.err('c_union_gap_nonneg_violation', max(0.0, PEh - min(1.0, sum_p)), EXACT_TOL, dump)
    # "G_e = 0 iff every part has at most one unit": strict when e has no p=1 unit; reported separately otherwise
    iff_ok = (abs(G) <= 1e-12) == all(len(u) <= 1 for _, u in parts)
    has_p1 = any(p[t] >= 1.0 for t in e)
    tr.flag('c_iff_G0_iff_parts_le_1_no_p1_unit' if not has_p1 else 'c_iff_G0_iff_parts_le_1_event_has_p1_unit',
            iff_ok, dict(dump, G=G), strict=not has_p1)


def check_E6(n_inst, rng):
    tr = Tracker()
    n_p1, n_ev = 0, 0
    for i in range(n_inst):
        N = int(rng.integers(4, 13))
        p = rng.uniform(0.02, 0.98, N)
        has_p1 = (i % 4 == 3)                                                # every 4th instance has p_t = 1 unit(s)
        if has_p1:
            p[rng.choice(N, size=int(rng.integers(1, 3)), replace=False)] = 1.0
            n_p1 += 1
        blocks = excl_partition(p, rng.random(N) < 0.2)
        inst = dict(N=N, p=p, blocks=blocks, mode='p1' if has_p1 else 'std', rng=rng)
        PE, PI = e6_laws_checks(tr, inst)
        inst = {k: v for k, v in inst.items() if k != 'rng'}
        for e in e6_events(rng, N, p, blocks):
            e6_event_checks(tr, inst, PE, PI, e, 'ev')
            n_ev += 1
    return tr.result(n_inst, n_events=n_ev, n_instances_with_p_eq_1=n_p1)


# --------------------------------------------------------------------------- E7
E7_INSTANCES = [
    # (name, p, forced start indices)
    ('N12', [0.30, 0.45, 0.25, 0.60, 0.15, 0.50, 0.35, 0.40, 0.20, 0.55, 0.10, 0.70], [0, 5]),
    ('N14_with_p1', [0.20, 0.35, 0.50, 0.90, 0.25, 0.30, 1.00, 0.40, 0.45, 0.15, 0.60, 0.55, 0.05, 0.80], [0, 3, 9]),
]


def z_score(freq, q, R):
    if q <= 0.0 or q >= 1.0:
        return 0.0 if abs(freq - q) < 1e-15 else float('inf')
    return (freq - q) / np.sqrt(q * (1 - q) / R)


def check_E7(seed, R):
    tr = Tracker()
    try:
        if PLUGIN_DIR not in sys.path:
            sys.path.insert(0, PLUGIN_DIR)
        from samplers import simulate_excl_openloop as simulate_excl_openloop
    except Exception as ex:                                                  # import problem is a failure, not a skip
        tr.flag('import_simulate_excl_openloop', False, dict(error=repr(ex)))
        return tr.result(0)
    zmax = {}
    for i, (name, pl, starts) in enumerate(E7_INSTANCES):
        p = np.array(pl)
        N = len(p)
        forced = np.zeros(N, bool)
        forced[starts] = True
        blocks = excl_partition(p, forced)                                   # independent re-computation of the blocks
        out = simulate_excl_openloop(p, forced, R, np.random.default_rng([seed, i]), period_len=6)
        S = out['S']
        new = np.zeros(N, bool)
        new[[lo for lo, _ in blocks]] = True
        tr.flag('blocks_agree_with_plugin_NEW', bool((out['NEW'] == new[None, :]).all()), dict(name=name, blocks=blocks))
        starts_idx = np.array([lo for lo, _ in blocks])
        tr.flag('at_most_one_upload_per_block', int(np.add.reduceat(S.astype(np.int32), starts_idx, axis=1).max()) <= 1,
                dict(name=name))
        zs = [z_score(S[:, t].mean(), p[t], R) for t in range(N)]
        ev_rng = np.random.default_rng([seed, 70 + i])
        events = [random_event(ev_rng, N, p, k) for k in ('contig', 'contig', 'contig', 'subset', 'subset', 'subset',
                                                          'heavy', 'heavy')]
        events += [w for w in (block_subset_event(ev_rng, blocks) for _ in range(2)) if w is not None]
        zev = []
        for e in events:
            _, ms = union_bound_Ue(blocks, p, e)
            q = 1.0 - float(np.prod([1 - m for m in ms]))
            zev.append(z_score(S[:, e].any(1).mean(), q, R))
        dump = dict(name=name, p=p, blocks=blocks, R=R, z_units=zs, events=events, z_events=zev)
        tr.err('max_abs_z_units', np.max(np.abs(zs)), 4.5, dump)
        tr.err('max_abs_z_events', np.max(np.abs(zev)), 4.5, dump)
        zmax[name] = dict(units=float(np.max(np.abs(zs))), events=float(np.max(np.abs(zev))), n_events=len(events))
    return tr.result(len(E7_INSTANCES), R=R, max_abs_z=zmax)


# --------------------------------------------------------------------------- main
CHECK_NAMES = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7']
DESCR = dict(E1='Prop. 1 (bound and attainment): LP max = U_e, EXCL attains, W=0', E2='Prop. 1 (decomposition): E[sumM]-P(hit)=R+W on LP vertices',
             E3='Prop. 1 (single start / delayed information): marginals-only LP = min(1,sum p); delayed info exceeds U_e',
             E4='Prop. 1 (conditional marginals): violating designs have a biased {0,1} predictor', E5='Prop. 1 (closed-loop counterexample): counterexample numbers',
             E6='Cor. 1: exact EXCL vs independent laws', E7='simulate_excl_openloop vs exact law (z-scores)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n-inst', type=int, default=200)
    ap.add_argument('--seed', type=int, default=20261009)
    ap.add_argument('--only', type=str, default='', help='comma list, e.g. E1,E3')
    ap.add_argument('--r7', type=int, default=200000, help='number of runs R in E7')
    ap.add_argument('--smoke', action='store_true', help='n-inst 5, E7 with R=20000, writes results/enum_smoke.json')
    args = ap.parse_args()
    if args.smoke:
        args.n_inst, args.r7 = 5, 20000
    names = [s.strip().upper() for s in args.only.split(',') if s.strip()] or CHECK_NAMES
    bad = [s for s in names if s not in CHECK_NAMES]
    if bad:
        ap.error('unknown check(s): %s' % bad)
    funcs = dict(E1=check_E1, E2=check_E2, E3=check_E3, E4=check_E4, E5=check_E5, E6=check_E6)
    results, rows = {}, []
    for k, name in enumerate(CHECK_NAMES):
        if name not in names:
            continue
        t0 = time.time()
        if name == 'E7':
            res = check_E7(args.seed, args.r7)
        else:
            # E4 re-generates E3's designs: same seed stream index as E3
            rng = np.random.default_rng([args.seed, CHECK_NAMES.index('E3') if name == 'E4' else k])
            res = funcs[name](args.n_inst, rng)
        res['seconds'] = time.time() - t0
        res['description'] = DESCR[name]
        results[name] = res
        rows.append((name, res))
        print('done %s in %.1fs' % (name, res['seconds']), file=sys.stderr)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, 'enum_smoke.json' if args.smoke else 'enum_results.json')
    meta = dict(n_inst=args.n_inst, seed=args.seed, r7=args.r7, smoke=args.smoke, only=names)
    with open(out_path, 'w') as fh:
        json.dump(jsonable(dict(meta=meta, results=results)), fh, indent=1, default=str)
    print('\n| check | n inst | max abs err | strict failures | pass | s | what |')
    print('|---|---|---|---|---|---|---|')
    for name, r in rows:
        print('| %s | %d | %.3g | %d | %s | %.1f | %s |' % (name, r['n_instances'], r['max_abs_err'],
                                                         r['n_strict_failures'], 'PASS' if r['pass_'] else 'FAIL',
                                                         r['seconds'], r['description']))
    print('\nJSON: %s' % os.path.relpath(out_path, HERE))
    return 0 if all(r['pass_'] for _, r in rows) else 1


if __name__ == '__main__':
    sys.exit(main())
