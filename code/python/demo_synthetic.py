"""Synthetic demonstration (no real data): independent draws, EXCL and Madow at the same Active Measurement probabilities.

A clustered population is generated (recordings, activity bouts of several consecutive units, a noisy score). Probabilities come from the
open-loop rule; the three designs differ only in the draw. Prints coverage of the nominal 95% interval, event recall, nominal and matched width
relative to independent draws, and the realised upload rate. Takes well under a minute.   python demo_synthetic.py [--R 400] [--b 0.05]
"""
import argparse, time
import numpy as np
import samplers as S
import estimator as E
import events as V


def make_population(n, rng, mean_rec=500, p_on=0.04, p_off=0.15, snr=1.3):
    """Recordings of random length; activity is a two-state Markov chain (bouts of about 1/p_off units); the score is the activity plus
    unit noise, smoothed over three units."""
    gaps = np.zeros(n, bool)
    t = 0
    while t < n:
        t += int(rng.integers(mean_rec // 2, 3 * mean_rec // 2))
        if t < n: gaps[t] = True
    y = np.zeros(n); on = 0
    for t in range(n):
        if gaps[t] or t == 0: on = 0
        on = (rng.random() < p_on) if not on else (rng.random() > p_off)
        y[t] = on
    raw = snr * y + rng.normal(size=n)
    s = np.convolve(raw, np.ones(3) / 3, mode='same')
    return s, y, gaps


def run(R=400, b=0.05, seed=1, n_test=8000, n_dev=8000):
    rng = np.random.default_rng(seed)
    s_dev, y_dev, _ = make_population(n_dev, rng)
    s, y, gaps = make_population(n_test, rng)
    cal = S.calibrate(s_dev, y_dev)
    jb = S.bin_of(cal, s); q = cal['q'][jb]
    q_hot = float(y_dev.mean())                                   # onset alignment level: development prevalence
    forced = S.block_flags(q, q_hot, 2, gaps)
    p = S.external_probs('AM', s, cal, b, s_dev)                  # open-loop probabilities (Active Measurement)
    events = V.events_from_units(y, gaps)
    theta = float(y.mean()); rec_unit = np.cumsum(np.where(np.arange(n_test) == 0, True, gaps)) - 1
    print('population: %d units, %d recordings, prevalence %.3f, %d events, mean probability %.4f (b = %.2f), %d replications'
          % (n_test, int(gaps.sum()) + 1, theta, len(events), p.mean(), b, R))
    out = {}
    t0 = time.time()
    # independent draws, fixed periods of 200 units
    Sd, Pd = S.draw_independent(p, R, np.random.default_rng(seed + 1))
    est, vh, vp, _ = E.prequential(y, Sd, Pd, jb, cal, 200, 0.3, 'intercept')
    out['Independent'] = (est, np.sqrt(np.maximum(np.maximum(vh, vp), 0)), Sd)
    # EXCL
    o = S.simulate_excl_openloop(p, forced, R, np.random.default_rng(seed + 2), 200)
    est, vh, vp = E.preq_blocks(y, o['S'], o['P'], jb, cal, o['PID'], 0.3)
    out['EXCL'] = (est, np.sqrt(np.maximum(np.maximum(vh, vp), 0)), o['S'])
    # Madow (single start per recording)
    o = S.simulate_madow(p, gaps, forced, R, np.random.default_rng(seed + 3), 200)
    est, vh, vp = E.preq_blocks(y, o['S'], o['P'], jb, cal, o['PID'], 0.3)
    out['Madow (single start)'] = (est, np.sqrt(np.maximum(np.maximum(vh, vp), 0)), o['S'])
    w_ref = E.cm_width(theta, *out['Independent'][:2])
    print('\n%-22s %9s %13s %14s %14s %9s' % ('design', 'coverage', 'event recall', 'nominal width', 'matched width', 'rate/b'))
    for name, (est, se, Sx) in out.items():
        cov = E.coverage(theta, est, se); rec = V.event_metrics(Sx, y, events)['event_recall']
        nominal = float(np.mean(2 * E.Z95 * se)) / float(np.mean(2 * E.Z95 * out['Independent'][1]))
        mw = E.cm_width(theta, est, se) / w_ref
        print('%-22s %9.3f %13.3f %14.3f %14.3f %9.3f' % (name, cov, rec, nominal, mw, Sx.mean() / b))
    ub = float(np.mean([min(1.0, p[a:c].sum()) for a, c in events]))
    print('\nunion bound on event recall at these probabilities: %.3f' % ub)
    print('elapsed %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--R', type=int, default=400); ap.add_argument('--b', type=float, default=0.05)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args(); run(a.R, a.b, a.seed)
