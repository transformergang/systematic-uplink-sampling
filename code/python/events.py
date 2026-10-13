"""Events (maximal runs of positive units, gaps of one negative unit merged, never across a recording start) and event recall."""
import numpy as np


def events_from_units(y, gaps=None, intervals=None, ev_thr=0.5, merge_gap=1):
    """(k,2) int array of [start,end) unit-index intervals. Runs of positive units merged across at most merge_gap negative units, never across a recording start (gaps[t]=True: unit t starts a new recording)."""
    if intervals is not None:
        return np.asarray(intervals, int).reshape(-1, 2)
    pos = np.asarray(y) >= ev_thr; n = len(pos)
    gaps = np.zeros(n, bool) if gaps is None else np.asarray(gaps, bool)
    runs, a = [], None
    for t in range(n):
        if a is not None and (not pos[t] or gaps[t]):
            runs.append([a, t]); a = None
        if pos[t] and a is None:
            a = t
    if a is not None: runs.append([a, n])
    ev = []
    for r in runs:
        if ev:
            gap_len = r[0] - ev[-1][1]
            if 0 < gap_len <= merge_gap and not gaps[ev[-1][1]:r[0] + 1].any():
                ev[-1][1] = r[1]; continue
        ev.append(list(r))
    return np.array(ev, int).reshape(-1, 2)


def event_metrics_runs(S, y, events, ev_thr=0.5):
    """Per-run arrays: window recall, event recall (fraction of events with at least one upload) and median delay of the first upload."""
    R = S.shape[0]; pos = np.asarray(y) >= ev_thr
    Sp = S & pos[None, :]
    win = S[:, pos].mean(1) if pos.any() else np.full(R, np.nan)
    k = len(events)
    if k == 0:
        return dict(win_recall=win, event_recall=np.full(R, np.nan), delay_med=np.full(R, np.nan))
    caught = np.zeros((R, k), bool); delay = np.full((R, k), np.nan)
    for i, (a, c) in enumerate(events):
        blk = Sp[:, a:c]; h = blk.any(1)
        caught[:, i] = h; delay[h, i] = np.argmax(blk[h], 1)
    with np.errstate(all='ignore'):
        dm = np.array([np.median(d[~np.isnan(d)]) if (~np.isnan(d)).any() else np.nan for d in delay])
    return dict(win_recall=win, event_recall=caught.mean(1), delay_med=dm)


def event_metrics(S, y, events, ev_thr=0.5):
    """Averages over the replications: event recall, miss rate, median delay, window recall."""
    m = event_metrics_runs(S, y, events, ev_thr)
    er = float(np.mean(m['event_recall']))
    return dict(win_recall=float(np.mean(m['win_recall'])), event_recall=er, event_miss=1 - er,
                delay_med=None if np.isnan(m['delay_med']).all() else float(np.nanmean(m['delay_med'])), n_events=int(len(events)))
